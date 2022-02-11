#######################################################################################
# Function to run OPENMALARIA                                                         #
#                                                                                     #
# Task:  run each simulations                                                         #
# Input: XML files                                                                    #
# Outputs: One Outupt files for each simulation                                       #
#                                                                                     #
# authors: thiery.masserey@swisstph.ch adapted from andrewjames.shattock@swisstph.ch  #
#######################################################################################


# load the function needed
require(plyr)
require(dplyr)
require(readr)
require(stringr)
require(tidyr)

# --------------------------------------------
# Run OpenMalaria for all generated xml files.
# --------------------------------------------
run_model = function(pm,
                     n_jobs,
                     xml_path = NA,
                     sim_path = NA) {
  # message in the console
  message("  - Running OpenMalaria")
  
  # Create a new log file for the cluster jobs
  log_file = file.path(pm$pth$log_files, "scicore_log.txt")
  create_bash_log(log_file)
  
  # Construct slurm array command for running in parallel
  slurm_array = paste0("--array=1-", n_jobs, "%250")
  
  # Add country sub directory to file paths
  xml_path = paste0(xml_path, file_sep())
  sim_path = paste0(sim_path, file_sep())
  
  # Concatenate system command
  sys_command = paste(
    "sbatch",
    slurm_array,
    "bash_openmalaria.sh",
    xml_path,
    sim_path,
    pm$pth$om_files,
    log_file
  )
  
  # Invoke this command
  system(sys_command)
  
  # Wait for all cluster jobs to complete
  wait_for_jobs(pm, log_file, n_jobs)
}