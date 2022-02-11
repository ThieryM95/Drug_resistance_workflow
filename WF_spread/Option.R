######################################################################################
# Option.R                                                                           #
#                                                                                    #
# Task: Set key simulation options: drug archetype (A, B or A+B), parameter spaces   #
# constrained parameter, number of simulations and seeds, directories                #
#                                                                                    #
# authors: thiery.masserey@swisstph.ch adapted from andrewjames.shattock@swisstph.ch #                                      #
######################################################################################


# download function in the environments
source("myRfunctions.R")
require(plyr)
require(dplyr)
require(xlsx)

# ------------------------
# Set options for analysis.
# -------------------------
set_options <- function(sample_num = 0, quiet = FALSE) {

  #-----general information------

  # Define analyze name
  pm <- list()
  pm$analysis_name <- "v2"

  # Define the user name (scicore login info)
  pm$user <- "masthi00"

  # define iteration (for adaptive sampling: the number update automatically at each round of adaptative sampling)
  pm$sample_num <- sample_num

  # ---- Sampling options ----

  # Number of latin hypercube samples to generate (prior to GP)
  opts <- list()
  opts$lhc_samples <- 250

  # Number of different seeds to simulate
  opts$n_seeds <- 3

  # Flag for resampling EIR values
  opts$do_resample <- TRUE

  # define drug archetype we want to model (will then choose the good xml, and parameter to vary, and parameter space)
  opts$drug <- "short" #  option: short (drug A), long (drug B), ACT_SR (drug A + B)

  # define if we want constrain some factors in the global sensitivity analysis
  opts$Type_analysis <- "Constrained" # option: Constrained (setting factors are not constrained, see constrained_parameter bellow), Not_constrained (setting factors are constrained see constrained_parameter bellow)

  # define the outcome of interest
  opts$om_outcome <- "Spread" # estimate the rate of spread (no other option)

  # ---- Gaussian process options ----

  # Select GP kernel function
  opts$gp_kernel <- "Matern5_2" #  other option: Gaussian, Matern5_2, Matern3_2 (Matern5_2 was selected as better fit)

  # Proportion of data withheld from training set to testing set
  opts$gp_test_split <- 0.2 #

  # Maximum number of iteration of optimisation algorithm
  opts$gp_max_iter <- 10000

  # ---- Adaptive sampling options ----

  # Maximum number of adaptive sampling attempts
  opts$sampling_max_iter <- 10

  # Number of points to re-simulate per arm (per setting)
  opts$n_adaptive_samples <- 250

  # Quantitative threshold for accepting GP performance
  opts$stop_criteria <- 0.975 # (not implemented yet: but it could be a threshold for the correlation coefficient betwen predicted and observe selection coefficient to  defin if we can stop adaptive sampling)

  # ---- Option for sensitivity analysis ----

  # number of sampling points
  opts$sa_n <- 150000

  # ---- finale the opts liste based on what was defined ----

  # Create options list
  pm$opts <- opts

  # Specify parameter values for constrained parameter
  pm <- constrained_parameter(pm, opts$Type_analysis, opts$drug) # see function bellow

  # Specify parameter range for not constrained parameter
  pm <- variable_parameter(pm, opts$drug, opts$Type_analysis) # see function bellow

  # define directory, and create folders
  pm <- set_dirs(pm, pm$opts$drug, pm$opts$Type_analysis) # see directories.R

  # ---- Display key options and return ----

  # Estimate the number of scenarios defined to inform the user
  n_param <- 0
  for (i in 1:length(pm$settings)) {
    n_param[i] <- length(pm$settings[[i]])
  }

  n_settings <- prod(n_param)
  pm$opts$n_total_jobs <- pm$opts$lhc_samples * n_settings * pm$opts$n_seeds

  # Only display is flag is on
  if (quiet == FALSE) {
    message(
      " - Total number of arm: ",
      format(n_settings, big.mark = ",")
    )

    message(
      " - Total number of simulations: ",
      format(pm$opts$n_total_jobs, big.mark = ",")
    )
  }
  return(pm)
}


# -----------------------------------------------------
# Define the parameter values of constrained parameters.
# ------------------------------------------------------
constrained_parameter <- function(pm, type_analysis, drug) {

  # creat a list vector
  settings <- list()

  # define the seasonality patern via montlhy EIR values
  sesonality1 <- c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
  sesonality2 <- c(5.9, 2, 8.8, 26.5, 63.8, 133.1, 230.5, 335.4, 411.3, 346.6, 189.1, 45.4)

  # define the value of constrained parameter in the case the analysis is contrained
  if (type_analysis == "Constrained") {

    # constrained parameter values when simulated short and long acting drug
    if (drug == "short" | drug == "long") {

      # Read in country data from input spreadsheet
      EIR <- c(5, 10, 500) # define eir simulated NB: do not go below 5 with seasonality pattern
      Seasonality <- c("sesonality1", "sesonality2") # define seasonality pattern simulated
      Access <- c(0.04, 0.5) # define probability to access treatment
      Resistance_Level <- c(7, 18) # define the level of resistance of the resistant gneotype to the short/long acting drug
      Dosage <- c(0, 1) # define treatment adherence (0: 60% adherence , 1; 100% adherence)

      # Merge all the information
      setting_data <- Reduce(merge, list(as.data.frame(EIR), as.data.frame(Access), as.data.frame(Resistance_Level), as.data.frame(Dosage), as.data.frame(Seasonality)))
      colnames(setting_data) <- c("eir", "Access", "Resistance_Level", "Dosage", "seasonality") #

      # Check that  there is no duplicate
      settings$eir <- unique(setting_data$eir)
      settings$Access <- unique(setting_data$Access)
      settings$Resistance_Level <- unique(setting_data$Resistance_Level)
      settings$Dosage <- unique(setting_data$Dosage)
      settings$seasonality <- data.frame(sesonality1, sesonality2)
    }
    

    # constrained parameter values when simulated ACT
    if (drug == "ACT_SR") {
      EIR <- c(5, 10, 500) # define eir simulated NB: do not go below 5 with seasonality pattern
      Seasonality <- c("sesonality1", "sesonality2") # define seasonality pattern simulated
      Access <- c(0.04, 0.5) # define probability to access treatment
      Resistance_Level <- c(7, 18) # define the level of resistance of the resistant genotype to the short acting drug
      Resistance_Level_long <- c(10) # define the level of resistance of the resistant genotype to the long acting drug
      Dosage <- c(0, 1) # define treatment adherence (0: 60% adherence , 1; 100% adherence)

      # Merge all the information
      setting_data <- Reduce(merge, list(as.data.frame(EIR), as.data.frame(Access), as.data.frame(Resistance_Level_long), as.data.frame(Resistance_Level), as.data.frame(Dosage), as.data.frame(Seasonality)))
      colnames(setting_data) <- c("eir", "Access", "Resistance_Level_long", "Resistance_Level", "Dosage", "seasonality") # remove eir before c("eir", "seasonality") ######################################

      # Check that  there is no duplicate
      settings$eir <- unique(setting_data$eir)
      settings$Access <- unique(setting_data$Access)
      settings$Resistance_Level_long <- unique(setting_data$Resistance_Level_long)
      settings$Resistance_Level <- unique(setting_data$Resistance_Level)
      settings$Dosage <- unique(setting_data$Dosage)
      settings$seasonality <- data.frame(sesonality1, sesonality2)
    }
  }

  # define the value of seasonality  in the case the analysis is notcontrained
  if (type_analysis == "Not_constrained") {
    Seasonality <- c("sesonality1", "sesonality2") # define the different seasonality pattern

    # Merge all the information
    setting_data <- data.frame(Seasonality)
    colnames(setting_data) <- c("seasonality")
    settings$seasonality <- data.frame(sesonality1, sesonality2)
  }

  # Append settings to pm list
  pm$settings <- settings

  return(pm)
}

# ------------------------------------------------------------------
# Define the parameter space for parameters that are not constrained
# ------------------------------------------------------------------
variable_parameter <- function(pm, drug, type_analysis) {

  # define parameters space in the case of constrained analysis
  if (type_analysis == "Constrained") {

    # define paramter space for short acting drug
    if (drug == "short") {

      # parameter name
      Parameter <- c("Fitness", "IC50_S", "half_life", "MKR", "Diagnostic")

      # maximum values
      max <- c(1, 0.009, 0.175, 31, 50)

      # minimum values
      min <- c(0.9, 0.0016, 0.035, 27.6, 2)
    }

    # define paramter space for long acting drug
    if (drug == "long") {

      # parameter name
      Parameter <- c("Fitness", "IC50_S", "half_life", "MKR", "Diagnostic", "Dosage_2")

      # maximum values
      max <- c(1, 0.03, 21, 5, 50, 40)

      # minimum values
      min <- c(0.9, 0.01, 7, 3.45, 2, 30)
    }
    

    # define paramter space for ACT
    if (drug == "ACT_SR") {

      # parameter names
      Parameter <- c("Fitness", "IC50_S_long", "half_life_long", "MKR_long", "Diagnostic", "Dosage_2_long", "IC50_S_short", "half_life_short", "MKR_short")

      # maximum values
      max <- c(1, 0.03, 21, 5, 50, 40, 0.009, 0.175, 31)

      # minimum values
      min <- c(0.9, 0.01, 7, 3.45, 2, 30, 0.0016, 0.035, 27.6)
    }
  }

  # define parameters space in the case of not_constrained analysis
  if (type_analysis == "Not_constrained") {

    # define paramter space for short acting drug
    if (drug == "short") {

      # parameter name
      Parameter <- c("Fitness", "Access", "eir", "IC50_S", "Resistance_Level", "half_life", "MKR", "Diagnostic")

      # maximum values
      max <- c(1, 0.55, 500, 0.009, 50, 0.175, 31, 50)

      # minimum values
      min <- c(0.9, 0.04, 5, 0.0016, 1, 0.035, 27.6, 2)
    }

    # define paramter space for long acting drug
    if (drug == "long") {

      # parameter names
      Parameter <- c("Fitness", "Access", "eir", "IC50_S", "Resistance_Level", "half_life", "MKR", "Diagnostic", "Dosage_2")

      # maximum values
      max <- c(1, 0.5, 500, 0.03, 18, 21, 5, 50, 40)

      # minimum values
      min <- c(0.9, 0.04, 5, 0.01, 1, 7, 3.45, 2, 30)
    }

    # define paramter space for ACT  drug
    if (drug == "ACT_SR") {

      # parameter names
      Parameter <- c("Fitness", "Access", "eir", "IC50_S_long", "Resistance_Level_long", "half_life_long", "MKR_long", "Diagnostic", "Dosage_2_long", "IC50_S_short", "Resistance_Level", "half_life_short", "MKR_short")

      # maximum values
      max <- c(1, 0.53, 500, 0.03, 18, 21, 5, 50, 40, 0.009, 18, 0.175, 31)

      # minimum values
      min <- c(0.9, 0.04, 5, 0.01, 1, 7, 3.45, 2, 30, 0.0016, 1, 0.035, 27.6)
    }
  }

  # creat a datafram of parameter and parameter space
  program_data <- data.frame(Parameter, max, min)

  # Convert dataframe to list
  prog <- as.list(program_data)

  # Names
  prog$prog_names <- Parameter

  # estimate the number of programs
  prog$n_progs <- length(prog$prog_names)


  # Append program details to pm list
  pm$prog <- prog

  return(pm)
}
