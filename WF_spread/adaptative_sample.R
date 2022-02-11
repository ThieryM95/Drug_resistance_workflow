######################################################################################
# ADAPTIVE SAMPLING                                                                  #
#                                                                                    #
# Functions to prepare for each adaptive sampling iteration.                         #
#                                                                                    #
# The Task:                                                                          #
# 1) Loop through the different arms to assess each GP model                         #
# 2) Check if adative sampling is necessary by assessing GP performance              #
# 3) Resample as necessary based on variance of GP emulator                          #
# 4) Generate new parameter tables and xml files                                     #
# 5) Rerun model and post-processing                                                 #
# 6) Repeat these steps until we're happy (or max iterations reached)                #
#                                                                                    #
# Input: Gaussian Process                                                            #
#                                                                                    #
# Output: Gaussian Process with better fit (improved with adaptative sampling)       #                                                #
#                                                                                    #
# authors: thiery.masserey@swisstph.ch adapted from andrewjames.shattock@swisstph.ch #                                      #
######################################################################################

# ----------------------------------------------------
# Parent function for all adaptive sampling processes.
# ----------------------------------------------------
run_adaptive_sampling <- function(pm) {

  # Creat a list vector to save the model and correlation of each round of GP for easy visualization
  outpout <- list()
  all_outpout <- list()

  # Repeat adaptive sampling UP TO sampling_max_iter times
  for (i in 1:(pm$opts$sampling_max_iter)) { # 1<-2
    
    # Define sample number
    pm$sample_num <- i
    
    # Message to console
    message("* Performing adaptive sampling iteration ", pm$sample_num)

    # ---- For each arm, if the accuracy of the gp is not good enough resample parameter  ----
    
    # Reset param table on each iteration
    as_param_table <- NULL

    # Loop through all of the different arms
    if (pm$opts$Type_analysis == "Not_constrained") {
      pm$settings$Dosage <- pm$settings$Access <- pm$settings$eir <- pm$settings$Resistance_Level <- 1
    }
    for (this_Dosage in pm$settings$Dosage) { 
      for (this_Access in pm$settings$Access) { 
        for (this_Resistance in pm$settings$Resistance_Level) { 
          for (this_eir in pm$settings$eir) { 
            for (this_season in names(pm$settings$seasonality)) { 

              # Define the parameter range that we can sample between
              cov_ranges <- matrix(c(pm$prog$min, pm$prog$max), ncol = 2)
              rownames(cov_ranges) <- pm$prog$prog_names

              # Construct file name based on arm values
              if (pm$opts$Type_analysis == "Constrained") {
                setting_name <- paste(pm$opts$om_outcome, this_season, "EIR", this_eir, "Access", this_Access, "Resistance", this_Resistance, "Dosage", this_Dosage, sep = "_") ### remove eir paste(pm$opts$om_outcome, "eir", this_eir, this_season, sep = "_")
                sample_name <- paste(setting_name, "sample", (pm$sample_num - 1), sep = "_")
              } else {
                setting_name <- paste(pm$opts$om_outcome, this_season, sep = "_") ### remove eir paste(pm$opts$om_outcome, "eir", this_eir, this_season, sep = "_")
                sample_name <- paste(setting_name, "sample", (pm$sample_num - 1), sep = "_")
              }

              # Load the GP of the arm
              lala <- paste0(pm$pth$gp_samples, sample_name, "_gp.RData")
              gp_file <- readRDS(lala)
              gp_model <- gp_file$gp_model

              # Save the content of gp
              outpout$setting_name <- gp_file
              names(outpout)[names(outpout) == "setting_name"] <- setting_name

              # Estimate the correlation between the predicted and true selection coefficient of the testing dataset
              CORELATION <- cor(gp_file$predict_test, gp_file$actual_test, use = "pairwise.complete.obs")

              # If correlation is low, have another round of adaptive sampling
              if (CORELATION <= pm$opts$stop_criteria) {

                # Sample the parameters in the region of the parameter space in which we are less confident in the prediction
                new_samples <- resample_points(pm, gp_model, cov_ranges, pm$opts$drug, pm$opts$Type_analysis, this_Resistance, this_Dosage) # see function bellow

                # Concatenate into new parameter table that contain simulation for all arms
                as_param_table <- bind_param_table(as_param_table, new_samples, this_season, this_eir, this_Access, this_Resistance, this_Dosage, pm$opts$Type_analysis) # see function bellow
                
              }
            }
          }
        }
      }
    }

    # Save the GP fit for visualizing the result in an easy way later
    all_outpout[[i]] <- outpout
    outpout <- list()

    # If Null table => all GP models satisfy stopping criteria => we're done, no more need to sample
    if (is.null(as_param_table)) break()
    
    # ---- Add finishing touches to parameter table ----

    # Redefine the name of the variable that are constrained
    colnames(as_param_table)[1] <- "seasonality"
    if (pm$opts$Type_analysis == "Constrained") {
      colnames(as_param_table)[2] <- "eir"
      colnames(as_param_table)[3] <- "Access"
      colnames(as_param_table)[4] <- "Resistance_Level"
      if (pm$opts$drug == "ACT_SR") {
        colnames(as_param_table)[5] <- "Resistance_Level_long"
        colnames(as_param_table)[6] <- "Dosage"
      } else {
        colnames(as_param_table)[5] <- "Dosage"
      }
    }

    # Append scenario name column to dataframe
    scenario_name <- paste0("scenario_", 1:nrow(as_param_table))
    as_param_table <- cbind(as.data.frame(scenario_name), as_param_table)


    # Repeat the scenarios for a predefined number of seeds
    as_param_table <- merge(as_param_table, data.frame(seed = 1:pm$opts$n_seeds))

    # Update number of jobs to run
    pm$opts$n_total_jobs <- nrow(as_param_table)

    # Message to the console
    message(
      " - Additional simulations: ",
      format(pm$opts$n_total_jobs, big.mark = ","))

    # ---- Redo all processes, including GP fitting ----

    # Create new sub-directories for latest sample number
    pm <- set_dirs(pm, pm$opts$drug, pm$opts$Type_analysis)

    # 1) Generate paramter table
    generate_param_table(pm, param_table = as_param_table)

    # 2) Generate seed paterns
    n_jobs_unique <- simulation_sed_patterns(pm)
    message("   ~ Total number of scenarios: ", thou_sep(n_jobs_unique), " (with seed)")

    # 3) Generate xml file
    generate_xml_files(pm,
      n_jobs = n_jobs_unique,
      file_path = pm$pth$xml_base,
      sed_path = pm$pth$sim_sed,
      xml_path = pm$pth$sim_xml)

    # 4) Run OpenMalaria
    run_model(pm, n_jobs_unique,
      xml_path = pm$pth$sim_xml,
      sim_path = pm$pth$sim_out)

    # 5) Postprocessing and summary
    Postprocess(pm)
    
    # 6) Estimate selection coefficient
    SummaryResults(pm)

    # 7) Perform Guassian Process on model outcome
    run_gp(pm)

  }

  # Save all gp model and prediction
  gp_file <- paste0(pm$pth$gp_models, "All_gp.RData")
  saveRDS(all_outpout, file = gp_file)

  # Return to the user all gp model and prediction
  return(all_outpout)
}

# ----------------------------------------------------------------------------
# Function to Resamples a set of points to be rerun based on largest variance.
# ----------------------------------------------------------------------------
resample_points <- function(pm, gp_model, cov_ranges, drug, Type_analysis, this_Resistance, this_Dosage) {

  # Resample new points across the parameter space using LHC samplig
  new_points <- lhs(10000, cov_ranges)

  # Predict the output with the GP model on the new points
  if (drug == "short") {
    pred_obj <- predict(x = new_points, object = gp_model)
    new_points <- dataFrame(new_points, colNames = rownames(cov_ranges))
  }

  # If a long-acting drug is used need first to calculate the C_max for the prediction
  if (drug == "long" & Type_analysis == "Constrained") {

    # Creat a new variable C-max
    colnames(new_points) <- pm$prog$Parameter
    C_max <- rep(0, 10000)
    new_points <- cbind(new_points, C_max)

    # Calculate C_Max based on treatment adherence and dosage
    if (this_Dosage == 0) {
      new_points[, 7] <- (-8.368621e-17) + new_points[, 6] * (4.048628e-03)
    }
    if (this_Dosage == 1) {
      new_points[, 7] <- (-1.004235e-16) + new_points[, 6] * (5.440722e-03)
    }

    C_max_IC50 <- new_points[, 7] / (new_points[, 2])
    new_points <- cbind(new_points, C_max_IC50)
    
    # Predict the output with the GP model on the new points
    pred_obj <- predict(x = new_points[, c(1, 3:5, 8)], object = gp_model)

    # Convert to dataframe and append variance at these points
    new_points <- dataFrame(new_points)
  }

  # Same than above but for not constrained analysis
  if (drug == "long" & Type_analysis == "Not_constrained") {

    # Creat a new variable C_max
    colnames(new_points) <- pm$prog$Parameter
    C_max <- rep(0, 10000)
    new_points <- cbind(new_points, C_max)

    # Estimate Cmax base on treatment dosage
    new_points[, 10] <- (-1.004235e-16) + new_points[, 9] * (5.440722e-03)
    C_max_IC50 <- new_points[, 10] / (new_points[, 4])
    new_points <- cbind(new_points, C_max_IC50)

    # Predict the output with the GP model on the new points
    pred_obj <- predict(x = new_points[, c(1:3, 5:8, 11)], object = gp_model)

    # Convert to dataframe and append variance at these points
    new_points <- dataFrame(new_points)
  }

  # Same than above but for ACTs
  if (drug == "ACT_SR" & Type_analysis == "Not_constrained") {

    # creat a new variable C_max
    colnames(new_points) <- pm$prog$Parameter
    C_max <- rep(0, 1000)
    new_points <- cbind(new_points, C_max)

    # estimate C_max based on treatment dosage
    new_points[, 14] <- (-1.004235e-16) + new_points[, 9] * (5.440722e-03)

    C_max_IC50 <- new_points[, 14] / (new_points[, 4])
    new_points <- cbind(new_points, C_max_IC50)

    # Predict the output with the GP model on the new points
    pred_obj <- predict(x = new_points[, c(1:3, 5:8, 10:13, 15)], object = gp_model)

    # Convert to dataframe and append variance at these points
    new_points <- dataFrame(new_points)
  }

  # Same but for ACT if constrained analysis
  if (drug == "ACT_SR" & Type_analysis == "Constrained") { ### check the number

    # creat a new variable called C_max
    colnames(new_points) <- pm$prog$Parameter
    C_max <- rep(0, 10000)
    new_points <- cbind(new_points, C_max)

    # estimate C_max based on treatment dosage and adherence
    if (this_Dosage == 0) {
      new_points[, 10] <- (-8.368621e-17) + new_points[, 6] * (4.048628e-03)
    }
    if (this_Dosage == 1) {
      new_points[, 10] <- (-1.004235e-16) + new_points[, 6] * (5.440722e-03)
    }

    C_max_IC50 <- new_points[, 10] / (new_points[, 2] * pm$settings$Resistance_Level_long)
    new_points <- cbind(new_points, C_max_IC50)

    # Predict the output with the GP model on the new points
    pred_obj <- predict(x = new_points[, c(1, 3:5, 7:9, 11)], object = gp_model)

    # Convert to dataframe and append variance at these points
    new_points <- dataFrame(new_points)
  }

  # Order points by decreasing variance of the prediction to find the parameter space that we have more unsecurity in our prediciton
  new_points <- new_points[order(pred_obj$sd2 + pred_obj$nugs, decreasing = TRUE), ]

  # Select set of points with largest posterior predictive variance
  new_samples <- new_points[1:pm$opts$n_adaptive_samples, ]

  # Re-order to match parameter table
  new_samples <- new_samples[, pm$prog$prog_names]

  return(new_samples)
}

# ----------------------------------------------------------
# Add the constrained variable (arm) to the parameter table.
# ----------------------------------------------------------
bind_param_table <- function(as_param_table, new_samples, this_season, this_eir, this_Access, this_Resistance, this_Dosage, Type_analysis) {

  # for each analysis add the constrained variable to the new parameter table
  if (Type_analysis == "Not_constrained") {
    updated_table <- cbind.data.frame(this_season, new_samples, row.names = NULL)
  } else {
    if (pm$opts$drug == "ACT_SR") {
      updated_table <- cbind.data.frame(this_season, this_eir, this_Access, this_Resistance, pm$settings$Resistance_Level_long, this_Dosage, new_samples, row.names = NULL)
    } else {
      updated_table <- cbind.data.frame(this_season, this_eir, this_Access, this_Resistance, this_Dosage, new_samples, row.names = NULL)
    }
  }

  # save all parameters into the updated table
  updated_table <- rbind(as_param_table, updated_table)

  # return the table with all parameter needed
  return(updated_table)
}

# ---------------------------------------------------------------------
# Define the parameter table name for each adaptive sampling iteration.
# ---------------------------------------------------------------------
param_set_name <- function(sample_num = NA, all_samples = FALSE) {

  # Default name for all samples
  if (all_samples == TRUE) {
    name <- "parameter_table_all_samples"
  } else {

    # Set in one place so can be more easily changed
    name <- paste0("parameter_table_sample_", sample_num)
  }

  return(name)
}