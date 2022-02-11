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
set_dirs <- function(pm, drug) {

  # Initiate file path lists
  pth <- out <- list()

  # sample number
  sample <- paste0("sample_", pm$sample_num)

  # ---- Code and resource locations ----

  # Base path to code repositories
  base_stm <- file.path("/scicore/home/penny/masthi00/")

  # Parent paths to all input files relating to this project
  pth$code <- file.path(base_stm, "wf_esthablishment") ### working directory
  input <- file.path(pth$code, "SIM_FOLDER") ### file with input in the working floder

  # Path to output form scicore (if simulaiton run or not)
  pth$log_files <- file.path(pth$code, "JOB_OUT")

  # ---- Input files ----

  # Path to base line XML depending on the drug modeled
  if (drug == "short") {
    pth$xml_base <- file.path(input, "Esthablishment_short.xml")
  }
  if (drug == "long") {
    pth$xml_base <- file.path(input, "Esthablishment_long.xml")
  }
  if (drug == "ACT") {
    pth$xml_base <- file.path(input, "Esthablishment_ACT.xml")
  }
  
  # Path to OpenMalaria resource files
  pth$om_files <- file.path(base_stm, "OM_schema40_1") ## Path to OM


  # ---- Output directories and files ----

  # Parent path to all output files relating to this project
  file_stem <- file.path(base_stm, "OUT_Esthablishment", "Short_access_0.5_final")

  # Path to parameter table
  out$param_table <- file.path(file_stem, "0_parameters")

  # Paths to full simulation files
  out$sim_files <- file.path(file_stem, "1_sim_files")
  out$sim_sed <- file.path(out$sim_files, "sed_files", sample)
  out$sim_xml <- file.path(out$sim_files, "xml_files", sample)

  # Paths to post-processed files
  sim_output <- file.path(file_stem, "2_sim_outputs")
  out$sim_out <- file.path(sim_output, "raw_output", sample)
  out$processed <- file.path(sim_output, "processed", sample)
  out$summarised <- file.path(sim_output, "summarised")

  # Make all output directories
  make_out_dirs(out)

  # Append paths to pm list
  pm <- append_dirs(pm, pth, out)

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
