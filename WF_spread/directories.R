######################################################################################
# SET DIRECTORIES                                                                    #
#                                                                                    #
# Set and get directories in one place and creat any directory that do not currently #
# exist.                                                                             #
#                                                                                    #
# Outputs: A list of relevant directories (within pm$pth) which can be referenced    #
# elsewhere                                                                          #
#                                                                                    #
# authors: thiery.masserey@swisstph.ch adapted from andrewjames.shattock@swisstph.ch #                                      #
######################################################################################

# load the functions
source("myRfunctions.R")

# --------------------------------------------
# Define paths for project inputs and outputs.
# --------------------------------------------
set_dirs <- function(pm, drug, type_analysis) {

  # Initiate file path lists
  pth <- out <- list()

  # Sample number
  sample <- paste0("sample_", pm$sample_num)

  # ---- Code and resource locations ----

  # Base path to code repositories
  base_stm <- file.path("/scicore/home/penny/masthi00") # need to be modify with user name on scicore

  # Parent paths to all input files relating to this project
  pth$code <- file.path(base_stm, "WF_spread") # Working directory
  input <- file.path(pth$code, "SIM_FOLDER") # Folder that contains the XML files

  # ---- Select the good base XML files based on analysis----

  # base file if analysis is constrained
  if (type_analysis == "Constrained") {

    # basefile if model short acting drug
    if (drug == "short") {
      pth$xml_base <- file.path(input, "BASE_FAST_CSA.xml")
    }

    # basefile if model long acting drug
    if (drug == "long") {
      pth$xml_base <- file.path(input, "BASE_LONG_CSA.xml")
    }
    

    # basefile if model ACT
    if (drug == "ACT_SR") {
      pth$xml_base <- file.path(input, "ACT_Short_acting_resistance.xml")
    }
  }

  # base file if analysis is not constrained
  if (type_analysis == "Not_constrained") {

    # basefile if model short acting drug
    if (drug == "short") {
      pth$xml_base <- file.path(input, "BASE_FAST_GSA.xml")
    }

    # basefile if model long acting drug
    if (drug == "long") {
      pth$xml_base <- file.path(input, "BASE_LONG_GSA.xml")
    }

    # basefile if model ACT
    if (drug == "ACT_SR") {
      pth$xml_base <- file.path(input, "ACT_Short_acting_resistance_GSA.xml")
    }
  }

  # path to output form scicore (if simulaiton run or not)
  pth$log_files <- file.path(pth$code, "JOB_OUT")

  # Path to OpenMalaria resource files (OM version 40.1)
  pth$om_files <- file.path(base_stm, "OM_schema40_1") 

  # ---- Output directories and files ----

  # Parent path to all output files relating to this project
  file_stem <- file.path(base_stm, "OUT_spread", "FINAL_GSA_FAST_Seasonality")

  # Path to parameter table
  out$param_table <- file.path(file_stem, "0_parameters")

  # Paths to  simulation folders
  out$sim_files <- file.path(file_stem, "1_sim_files")

  # Paths to seed files
  out$sim_sed <- file.path(out$sim_files, "sed_files", sample)

  # Paths to xml files
  out$sim_xml <- file.path(out$sim_files, "xml_files", sample)

  # Paths to Outputs folders
  sim_output <- file.path(file_stem, "2_sim_outputs")

  # Path to OpenMalaria raw outputs
  out$sim_out <- file.path(sim_output, "raw_output", sample)

  # Path to Post processed Outputs
  out$processed <- file.path(sim_output, "processed", sample)

  # Path to summary outputs
  out$summarised <- file.path(sim_output, "summarised")

  # Paths to  GP fitted during every round of adaptative sampling
  out$gp_samples <- file.path(file_stem, "3_gp_models")

  # path to the GP fitted during the last round of adaptative sampling (use for the sensitivity analysis)
  out$gp_models <- file.path(file_stem, "3_gp_models", "models")

  # Paths to sensitivity analysis results
  out$results <- file.path(file_stem, "4_results")
  out$results_gp <- file.path(file_stem, "4_results")
  out$results_optim <- file.path(file_stem, "4_results")

  # Make all output directories
  make_out_dirs(out) # see function bellow

  # Append paths to pm list
  pm <- append_dirs(pm, pth, out) # see function bellow

  return(pm)
}


# ---------------------------------------------------------
# Make all output directories if they do not already exist.
# ---------------------------------------------------------
make_out_dirs <- function(out) {

  # Extract all path names in list
  pth_names <- names(out)

  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth <- out[[pth_name]]

    # If it does not already exist, create it
    if (!dir.exists(this_pth)) {
      dir.create(this_pth, recursive = TRUE)
    }
  } # Close directory loop
}

# ---------------------------------------------------------
# Concatenate separators and append directories to pm list.
# ---------------------------------------------------------
append_dirs <- function(pm, pth, out) {

  # Extract all path names in list
  pth_names <- names(out)

  # Loop through these path names
  for (pth_name in pth_names) {
    this_pth <- out[[pth_name]]

    # Add a file separator to end of output paths
    out[[pth_name]] <- paste0(this_pth, file_sep())
  }

  # Concatenate lists
  pm$pth <- c(pth, out)

  return(pm)
}
