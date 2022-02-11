######################################################################################
# Option.R                                                                           #
#                                                                                    #
# Task: Set key simulation options: drug archetype (A, B, A+B), number of simulation,#
# number of seeds, directories, constrain parameter                                  #
#                                                                                    #
# authors: thiery.masserey@swisstph.ch adapted from andrewjames.shattock@swisstph.ch #                                      #
######################################################################################


# Download function in the environments
source("myRfunctions.R")
require(plyr)
require(dplyr)
require(xlsx)

# ---------------------------------------------------------
# Set options for a locally called analysis.
# ---------------------------------------------------------
set_options <- function(do_step = NA, sample_num = 0, quiet = FALSE) {

  #-------general information------

  pm <- list()

  # Define the user name
  pm$user <- "masthi00"

  # Define analysis name
  pm$analysis_name <- "Esthablishment"

  # ---- Sampling options ----

  opts <- list()

  # Select the drug archetype we want to model
  opts$drug <- "short" #  option are long/short/ACT

  # Number of simulation (different mutation)
  opts$lhc_samples <- 10

  # Number of different seeds to simulate for each simulation 
  opts$n_seeds <- 30

  # ---- Post processing settings ----

  # Define the outcome of the analysis (only one option now)
  opts$om_outcome <- "esthablishment"

  # Define the sample iteration
  pm$sample_num <- sample_num # Adaptive sampling iteration number (not used in the workflow yet)

  # Create output directories
  pm$opts <- opts

  # Country specific settings and interventions
  pm <- constrained_parameter(pm) # see function bellow
  pm <- variable_parameter(pm, pm$opts$drug) # see function bellow

  # Define the directory
  pm <- set_dirs(pm, pm$opts$drug) # see directories.R

  # ---- Display key options and return ----

  # Number of scenarios defined - just to inform the user
  n_settings <- length(pm$settings$eir) * ncol(pm$settings$seasonality) * length(pm$settings$Access) * length(pm$settings$Dosage)
  pm$opts$n_total_jobs <- pm$opts$lhc_samples * n_settings * pm$opts$n_seeds

  # Only display is flag is on
  if (quiet == FALSE) {
    message(" - Running country: ", pm$country)
    message(" - Analysis name: ", pm$analysis_name)

    message(
      " - Total number of settings: ",
      format(n_settings, big.mark = ",")
    )

    message(
      " - Total number of simulations: ",
      format(pm$opts$n_total_jobs, big.mark = ",")
    )
  }

  return(pm)
}

# ---------------------------------------------------------
# Define country-specific setting details.
# ---------------------------------------------------------
constrained_parameter <- function(pm) {

  # Defined the constrained parameter values
  EIR <- c(5, 10, 500)
  Seasonality <- c("sesonality1")
  Access <- c(0.5)
  Dosage <- c(1)
  setting_data <- Reduce(merge, list(as.data.frame(EIR), as.data.frame(Access), as.data.frame(Dosage), as.data.frame(Seasonality)))

  # Merge the setting information
  colnames(setting_data) <- c("eir", "Access", "Dosage", "seasonality") # remove eir before c("eir", "seasonality") ######################################

  # Save the information into list setting
  settings <- list()
  settings$eir <- unique(setting_data$eir)
  settings$Access <- unique(setting_data$Access)
  settings$Dosage <- unique(setting_data$Dosage)

  # Seasonality profiles can be read in directly
  sesonality1 <- c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
  settings$seasonality <- data.frame(sesonality1)

  # Append settings to pm list
  pm$settings <- settings

  return(pm)
}

# ---------------------------------------------------------
# Set country-specific intervention details.
# ---------------------------------------------------------
variable_parameter <- function(pm, drug) {

  # Define the parameter investigate based on the drug archetype
  if (drug == "short") {
    Parameter <- c("Fitness", "Resistance_level", "IC50_S", "half_life", "MKR", "Diagnostic", "Importation")
    }

  if (drug == "long") {
    Parameter <- c("Fitness", "Resistance_level", "IC50_S", "half_life", "MKR", "Diagnostic", "Dosage_2", "Importation")
    }

  if (drug == "ACT") {
    Parameter <- c("Fitness", "Resistance_level", "Resistance_level_long", "IC50_S_long", "half_life_long", "MKR_long", "Diagnostic", "Dosage_2_long", "IC50_S_short", "half_life_short", "MKR_short")
  }

  # Merge information in one dataset
  program_data <- data.frame(Parameter)

  # Convert data frame to list
  prog <- as.list(program_data)

  # Names
  prog$prog_names <- Parameter

  # Easy access number of programs
  prog$n_progs <- length(prog$prog_names)

  # Append program details to pm list
  pm$prog <- prog

  return(pm)
}
