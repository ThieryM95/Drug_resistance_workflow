######################################################################################
# GAUSSIAN PROCESS                                                                   #
#                                                                                    #
# Task: Fit a GP to interpolate whole parameter space using OM simulations.          #
# Input: Summary table of parameter values and selection coefficient                 #
# Ouptut: GP fitted to each constrained arms                                         #
#                                                                                    #
# authors: thiery.masserey@swisstph.ch adapted from andrewjames.shattock@swisstph.ch #                                      #
######################################################################################

# Load function
require(hetGP)

# ----------------------------------------
# Parent function to perform all GP steps.
# ----------------------------------------
run_gp <- function(pm) {

  # Creat a list of elements
  indices <- list()

  # Message in the consoles
  message("* Fitting GP models")

  # Load the full summary data (summary table of each iteration)
  param_table_name <- param_set_name(all_samples = TRUE)
  summary_file <- paste0(param_table_name, "_", ".txt")
  summary_path <- paste0(pm$pth$summarised, summary_file)
  indicator_table <- read.table(summary_path, sep = "\t", header = TRUE, as.is = TRUE)

  # Calculate the ratio between maximum concentration and IC50 for long acting drug
  if (pm$opts$drug == "long" & pm$opts$Type_analysis == "Not_constrained") {
    indicator_table$C_max <- 0
    indicator_table$C_max <- (-1.004235e-16) + indicator_table$Dosage * (5.440722e-03) # 100% treatment adherence
    indicator_table$C_max_IC50 <- indicator_table$C_max / (indicator_table$IC50_S)
  }

  # calculate the ration between maximum concentration and IC50 for long acting drug
  if (pm$opts$drug == "long" & pm$opts$Type_analysis == "Constrained") {
    indicator_table$C_max <- 0
    indicator_table$C_max[indicator_table$Dosage == 0] <- (-8.368621e-17) + indicator_table$Dosage_2[indicator_table$Dosage == 0] * (4.048628e-03) # when 60% treatment adherence
    indicator_table$C_max[indicator_table$Dosage == 1] <- (-1.004235e-16) + indicator_table$Dosage_2[indicator_table$Dosage == 1] * (5.440722e-03) # when 100% treatment adherence
    indicator_table$C_max_IC50 <- indicator_table$C_max / (indicator_table$IC50_S)

    # need creat this parameter and define to one that the code can work for ACT later
    pm$settings$Resistance_Level_long <- 1
  }

  # calculate the ration between Max concentration and IC50 for long acting drug
  if (pm$opts$drug == "ACT_SR" & pm$opts$Type_analysis == "Not_constrained") {
    indicator_table$C_max <- 0
    indicator_table$C_max <- (-1.004235e-16) + indicator_table$Dosage_2_long * (5.440722e-03) # 100% treatment adherence
    indicator_table$C_max_IC50 <- indicator_table$C_max / (indicator_table$IC50_S_long)
  }

  # Calculate the ratio between Max concentration and IC50 for long acting drug
  if (pm$opts$drug == "ACT_SR" & pm$opts$Type_analysis == "Constrained") {
    indicator_table$C_max <- 0
    indicator_table$C_max[indicator_table$Dosage == 0] <- (-8.368621e-17) + indicator_table$Dosage_2_long[indicator_table$Dosage == 0] * (4.048628e-03)
    indicator_table$C_max[indicator_table$Dosage == 1] <- (-1.004235e-16) + indicator_table$Dosage_2_long[indicator_table$Dosage == 1] * (5.440722e-03)
    indicator_table$C_max_IC50 <- indicator_table$C_max / (indicator_table$IC50_S_long)
  }

  # Need creat this parameter and define to one that the code can work for ACT later
  if (pm$opts$drug == "short" & pm$opts$Type_analysis == "Constrained") {
    pm$settings$Resistance_Level_long <- 1
  }

  # Need creat these arms parameters and define the one that the code can work for ACT later
  if (pm$opts$Type_analysis == "Not_constrained") {
    pm$settings$Dosage <- pm$settings$Access <- pm$settings$eir <- pm$settings$Resistance_Level <- pm$settings$Resistance_Level_long <- 1
  }

  # Loop through all of the different arms
  for (this_Dosage in pm$settings$Dosage) { # 
    for (this_Access in pm$settings$Access) { #
      for (this_Resistance in pm$settings$Resistance_Level) { #
        for (this_eir in pm$settings$eir) { #
          for (this_season in names(pm$settings$seasonality)) { #
            for (this_Resistance_S in pm$settings$Resistance_Level_long) { #

              # ---- select the data specific to each arm -----
              
              # Construct file name based on arm values when constrained analysis
              if (pm$opts$Type_analysis == "Constrained") {
                setting_name <- paste(pm$opts$om_outcome, this_season, "EIR", this_eir, "Access", this_Access, "Resistance", this_Resistance, "Dosage", this_Dosage, sep = "_") ### remove eir paste(pm$opts$om_outcome, "eir", this_eir, this_season, sep = "_")
                sample_name <- paste(setting_name, "sample", pm$sample_num, sep = "_")

                # Reduce the table to only consider this arm
                setting_idx <- indicator_table$season == this_season &
                  indicator_table$Access == this_Access &
                  indicator_table$eir == this_eir &
                  indicator_table$Resistance_Level == this_Resistance &
                  indicator_table$Dosage == this_Dosage

                if (pm$opts$drug == "ACT_SR") {
                  setting_idx <- indicator_table$season == this_season &
                    indicator_table$Access == this_Access &
                    indicator_table$eir == this_eir &
                    indicator_table$Resistance_Level == this_Resistance &
                    indicator_table$Dosage == this_Dosage
                }
              
              # Construct file name based on arm values when not constrained analysis
              } else {
                setting_name <- paste(pm$opts$om_outcome, this_season, sep = "_") ### remove eir paste(pm$opts$om_outcome, "eir", this_eir, this_season, sep = "_")
                sample_name <- paste(setting_name, "sample", pm$sample_num, sep = "_")

                # Reduce the table to only consider this arm
                setting_idx <- indicator_table$season == this_season
              }

              # Message to console to know which arm is done
              message("  - ", setting_name)

              # select the scenario that fitt the arms
              setting_table <- indicator_table[setting_idx, ]

              # Remove annoying row numbers
              rownames(setting_table) <- NULL

              # Save all the scenario run in this arms
              setting_file <- paste0(pm$pth$gp_samples, sample_name, ".txt")
              write.table(setting_table, setting_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE) # save the table

              # ---- Train the model ----

              # Select input variables of interest
              all_progs <- pm$prog$prog_names ## select only paramter and indicators
              input_data <- setting_table

              # Train a GP model on this data using a cross validation approach
              gp_result <- split_train(pm, input_data, 0.2, pm$opts$drug, pm$opts$Type_analysis) # see function bellow

              # Save the fitted GP model to file - adaptive sample iteration
              gp_file <- paste0(pm$pth$gp_samples, sample_name, "Matern_3_gp.RData")
              saveRDS(gp_result, file = gp_file)

              # Save the fitted GP model to file - latest fitted model
              gp_file <- paste0(pm$pth$gp_models, setting_name, "_gp.RData")
              saveRDS(gp_result, file = gp_file)

              # save the fitted GP in a list to return to the user
              indices$sample_name <- gp_result
              names(indices)[names(indices) == "sample_name"] <- setting_name
            }
          }
        }
      }
    }
  }
  
  # return GP model of each arm 
  return(indices)
}

# ------------------------------------------------
# Simple split of train-test data to train the GP.
# ------------------------------------------------
split_train <- function(pm, input_data, p, drug, type_analysis) {

  # Define the parameter name
  all_progs <- pm$prog$prog_names

  # Modify parameter name (merge IC50 and Cmax in one parameter for long acting drug)
  if (drug == "long" & type_analysis == "Constrained") {
    all_progs <- c("Fitness", "half_life", "MKR", "Diagnostic", "C_max_IC50")
  }
  if (drug == "long" & type_analysis == "Not_constrained") {
    all_progs <- c("Fitness", "Access", "eir", "Resistance_Level", "half_life", "MKR", "Diagnostic", "C_max_IC50")
  }
  if (drug == "ACT_SR" & type_analysis == "Not_constrained") {
    all_progs <- c("Fitness", "Access", "eir", "Resistance_Level_long", "half_life_long", "MKR_long", "Diagnostic", "IC50_S_short", "Resistance_Level_short", "half_life_short", "MKR_short", "C_max_IC50")
  }
  if (drug == "ACT_SR" & type_analysis == "Constrained") {
    all_progs <- c("Fitness", "half_life_long", "MKR_long", "Diagnostic", "IC50_S_short", "half_life_short", "MKR_short", "C_max_IC50")
  }

  # Select all scenario with 1 seed
  input_data_1 <- input_data[input_data$seed == 1, ]

  # calculate the total number of scenario
  n_samples <- nrow(input_data_1)

  # Randomly sample the senario to be part of the test dataset
  n_test <- round(n_samples * p)
  test_idx <- sort(sample(1:n_samples, n_test))

  # Select the other scenairos as part of the training dataset
  train_idx <- setdiff(1:n_samples, test_idx)

  # Use the indices to seperate the data
  train_data <- input_data_1[train_idx, ]
  test_data <- input_data_1[test_idx, ]

  # Select all the other seed to be part of the good dataset
  train_data <- input_data[input_data$Fitness %in% train_data$Fitness, ]
  TEST_DATA <- input_data[input_data$Fitness %in% test_data$Fitness, ]

  # Select the  parameter columns + indicators columns only
  train_data <- train_data[, c(all_progs, "Indicator")]
  TEST_DATA <- TEST_DATA[, c(all_progs, "Indicator")]

  # Remove NA in the trained dataset
  train_data <- train_data[rowSums(is.na(train_data)) < 1, ]

  # For the test dataset calculate the mean selection coefficient (indicator) across seed
  test_data_2 <- test_data[FALSE, c(all_progs, "Indicator")]
  for (i in 1:length(test_data$Fitness)) {
    test_data_2[i, ] <- colMeans(TEST_DATA[test_data$Fitness[i] == TEST_DATA$Fitness, ], na.rm = TRUE)
  }

  # Train model using these training data
  trained_model <- gp_train(train_data) # see function bellow

  # Assess the model with both testing and training data
  gp_output <- gp_test(trained_model, train_data, test_data_2) # see function bellow

  # Append  to this output
  gp_output$gp_model <- trained_model
  
  # Return the gp model and prediction
  return(gp_output)
}

# -----------------------------------------------------------
# Train a GP regression with hetGP given a training data set.
# -----------------------------------------------------------
gp_train <- function(train_data) {

  # Define indices of parameters and response
  n_params <- ncol(train_data) - 1
  param_idx <- 1:n_params
  response_idx <- n_params + 1

  # Find and remove duplicates
  prep_data <- find_reps(
    X = as.matrix(train_data[, param_idx]),
    Z = as.matrix(train_data[, response_idx]),
    rescale = FALSE, normalize = FALSE)

  # Prepare input for mleHetGP model
  X <- list(
    X0 = as.matrix(prep_data$X0),
    Z0 = as.matrix(prep_data$Z0),
    mult = prep_data$mult)

  # Run GP model
  gp_model <- mleHetGP(
    X = X, Z = prep_data$Z,
    covtype = "Matern5_2")
  
  # Return gp model
  return(gp_model)
}

# -----------------------------------------------
# Run GP with test data to determine model error.
# -----------------------------------------------
gp_test <- function(gp_model, train_data, test_data_2) {

  # Define indices of parameters and response
  n_params <- ncol(train_data) - 1
  param_idx <- 1:n_params
  response_idx <- n_params + 1

  # Apply the model on the test data
  predict_test <- predict(x = as.matrix(test_data_2[, param_idx]), object = gp_model)
  predict_test <- predict_test$mean #### look at the mean of predictions

  # The actual response associated with the test data
  actual_test <- test_data_2[, response_idx]

  # Repeat this for the training data
  predict_train <- predict(x = as.matrix(train_data[, param_idx]), object = gp_model)
  predict_train <- predict_train$mean

  # The actual response associated with the test data
  actual_train <- train_data[, response_idx]

  # Concatenate output into list
  gp_output <- list(
    actual_test = actual_test,
    predict_test = predict_test,
    actual_train = actual_train,
    predict_train = predict_train)
  
  # Return the predicted data and true value
  return(gp_output)
}
