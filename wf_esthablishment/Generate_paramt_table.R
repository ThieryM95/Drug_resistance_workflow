######################################################################################
#      Generate the parameter tables for the simulation that will be run             #
#                                                                                    #
#                                                                                    #
# Task: This function create the table of input parameter by selecting               #
#       combination of parameters with a know selection coefficient. The parameter   #
#       combination are selected to cover the range of selection coefficient estimate#
#       in this setting                                                              #
#                                                                                    #
# Input: number of simulation, number of seed, constrained parameter values,         #
#        and table of input parameter combination and their selection coefficient    #
#                                                                                    #
# Output: A table of the parameter value of each parameter for each simulation       #
#                                                                                    #
# authors: thiery.masserey@swisstph.ch                                               #
######################################################################################


# --------------------------------------------
# Generate the input table for the simulation.
# --------------------------------------------
generate_param_table <- function(pm, param_table = NULL) {
  
  # Message in the Console
  message("* Generating parameter table")

  #---- Select combination of parameter to simulate ----

  # File name for this parameter table
  param_table_name <- param_set_name(sample_num = 0)
  param_file <- paste0(pm$pth$param_table, param_table_name, ".txt")

  # For each drug archetype, select the table containing the estimation of selection coefficient for different combination of parameter
  if (pm$opts$drug == "long") {
    Param <- read.table("/scicore/home/penny/masthi00//wf_esthablishment/SIM_FOLDER/parameter_table_all_samples_long.txt", header = T)
  }
  if (pm$opts$drug == "short") {
    Param <- read.table("/scicore/home/penny/masthi00//wf_esthablishment/SIM_FOLDER/parameter_table_all_samples_short.txt", header = T)
  }
  if (pm$opts$drug == "ACT") {
    Param <- read.table("/scicore/home/penny/masthi00//wf_esthablishment/SIM_FOLDER/parameter_table_all_samples_ACT.txt", header = T)
    colnames(Param)[4] <- "Resistance_Level"
  }

  PARAM <- list()

  # Do a loop across each setting
  for (this_Dosage in pm$settings$Dosage) {
    for (this_Access in pm$settings$Access) {
      for (this_eir in pm$settings$eir) {
        for (this_season in names(pm$settings$seasonality)) {


          # Select each simulation for which we have estimate of the selection coefficient in this setting
          Settings <- Param[Param$eir == this_eir & Param$Access == this_Access & Param$Dosage == this_Dosage & Param$seasonality == this_season, ]
          Settings_1seed <- Settings[Settings$seed == 1, ]

          # select all parameter depending on drug archetype
          all_progs <- pm$prog$prog_names

          if (pm$opts$drug == "long") {
            mean_settings <- Settings_1seed[, c(1:10, 12)]
          }
          if (pm$opts$drug == "short") {
            mean_settings <- Settings_1seed[, c(1:9, 11)]
          }

          if (pm$opts$drug == "ACT") {
            mean_settings <- Settings_1seed[, c(1:14, 16)]
          }

          # Calculate the mean selection coefficient across seed
          mean_settings$Indicator <- 0
          RL <- levels(as.factor(Settings_1seed$Resistance_Level))
          RL <- as.numeric(RL)
          for (this_RL in RL) {
            for (i in 1:length(Settings$Fitness)) {
              mean_settings$Indicator[mean_settings$Resistance_Level == this_RL & Settings_1seed$Fitness[i] == mean_settings$Fitnes] <- mean(Settings$Indicator[Settings_1seed$Fitness[i] == Settings$Fitness & Settings$Resistance_Level == this_RL], na.rm = T)
            }
          }

          # Look at the maximum and minimum selection coefficient in this setting
          mean_settings <- mean_settings[complete.cases(mean_settings), ]
          MAX <- max(mean_settings$Indicator)
          MIN <- min(mean_settings$Indicator)

          # Define the mean to be equal to zero (as do not want estimate the probability of establishment of mutation that are not selected)
          if (MIN <= 0) {
            MIN <- 0
          }

          # Define the sequence of selection coefficient investigate to cover the whole selection coefficient range
          Selected_IND <- seq(MIN, MAX, length.out = pm$opts$lhc_samples)

          # Select the combination of parameters that have the selection coefficient of interest
          Param_table <- mean_settings[1, ]
          for (i in 1:pm$opts$lhc_samples) {
            min_2 <- min((mean_settings$Indicator - Selected_IND[i])^2)
            Param_table[i, ] <- mean_settings[mean_settings$Indicator == (Selected_IND[i] - sqrt(min_2)) | mean_settings$Indicator == (sqrt(min_2) + Selected_IND[i]), ]
          }

          # save the simulation selected for this settings
          setting_name <- paste(this_season, "EIR", this_eir, "Access", this_Access, "Resistance", this_Dosage, sep = "_") ### remove eir paste(pm$opts$om_outcome, "eir", this_eir, this_season, sep = "_")
          PARAM[[setting_name]] <- Param_table
        }
      }
    }
  }


  # Estimated the number of setting
  n_settings <- length(pm$settings$eir) * ncol(pm$settings$seasonality) * length(pm$settings$Access) * length(pm$settings$Dosage)

  # Remove the list into a dataframe
  if (n_settings == 1) {
    PARAM_2 <- PARAM[[1]]
  } else {
    if (n_settings == 2) {
      PARAM_2 <- rbind(PARAM[[1]], PARAM[[2]])
    } else {
      PARAM_2 <- rbind(PARAM[[1]], PARAM[[2]])

      for (i in 3:n_settings) {
        PARAM_2 <- rbind(PARAM_2, PARAM[[i]])
      }
    }
  }

  #---- Estimate the importation rate in each setting to model a specific mutation rate ----

  # Download the data about the prevalence in each setting depending on the drug archetype
  if (pm$opts$drug == "short") {
    Prevalance <- read.table("/scicore/home/penny/masthi00/wf_esthablishment/SIM_FOLDER/Prevalance_2_short.txt", header = T)
  }
  if (pm$opts$drug == "long") {
    Prevalance <- read.table("/scicore/home/penny/masthi00/wf_esthablishment/SIM_FOLDER/Prevalance_2_long.txt", header = T)
  }
  if (pm$opts$drug == "ACT") {
    Prevalance <- read.table("/scicore/home/penny/masthi00/wf_esthablishment/SIM_FOLDER/Prevalance_2_ACT.txt", header = T)
  }


  # add a variable importation rate in the input table
  PARAM_2$Importation <- 0
  
  # define the number of generation per year
  g <- 6
  
  # define the mutation rate
  u <- 2 * 10^-6

  # Do a loop across each setting
  for (this_Dosage in pm$settings$Dosage) { # this_Dosage = pm$settings$Dosage[1]
    for (this_Access in pm$settings$Access) { # this_Access = pm$settings$Access[1]
      for (this_eir in pm$settings$eir) { # this_eir = pm$settings$eir[2]
        for (this_season in names(pm$settings$seasonality)) { # this_season = names(pm$settings$seasonality)[1]
          
          #  Estimate the impotation rate to model a specific mutaiton rate in this setting
          PARAM_2$Importation[PARAM_2$eir == this_eir & PARAM_2$Access == this_Access & PARAM_2$Dosage == this_Dosage & PARAM_2$seasonality == this_season] <- Prevalance$prevalance[Prevalance$eir == this_eir & Prevalance$Access == this_Access & Prevalance$Dosage == this_Dosage & Prevalance$seasonality == this_season] * u * g * 1000 * 2
        
          }
      }
    }
  }

  # add scenario name
  scenario_name <- paste("scenario", 1:nrow(PARAM_2), sep = "_")
  param_table <- cbind(scenario_name, PARAM_2)

  # add seed
  seed <- 1:pm$opts$n_seeds
  param_table <- merge(param_table, as.data.frame(seed))

  # Reorder the order of variable to make it more easy to generate the seed pattern
  if (pm$opts$drug == "long") {
    param_table_2 <- param_table[, c(1, 12, 2:11, 14, 15, 13)]
  }
  if (pm$opts$drug == "short") {
    param_table_2 <- param_table[, c(1, 11, 2:10, 13, 14, 12)]
  }

  if (pm$opts$drug == "ACT") {
    param_table_2 <- param_table[, c(1, 16, 2:15, 18, 19, 17)]
  }

  # Save the paramter table
  write.table(param_table_2,
    file = param_file, sep = "\t",
    quote = FALSE, col.names = TRUE, row.names = FALSE
  )
}
