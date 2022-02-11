#########################################################################################
# SIMULATION SETUP                                                                     ##
#                                                                                      ##
# TASK: For each simulations, creat the seed patern and xml file                       ##
# INPUTS: Parameter table                                                              ##
# OUTPUTS: Seed files and then XML file                                                ##
# Written by thiery.masserey@swisstph.ch adapted from andrewjames.shattock@swisstph.ch ##
#########################################################################################

# load funcitons requiered
source("myRfunctions.R")
require(pracma)
require(stringr)
require(tgp)
require(qdap)
require(xlsx)

# -------------------------------------------------------------------
# Generate xml seed replacement patterns for full country simulations.
# -------------------------------------------------------------------
simulation_sed_patterns <- function(pm) {

  # Message in the console
  message("  - Generating scenarios")

  # File path to  parameter table
  param_table_name <- param_set_name(sample_num = pm$sample_num)
  param_file <- paste0(pm$pth$param_table, param_table_name, ".txt")

  # Read the parameter table
  if (file.exists(param_file)) {
     param_table <- read.table(param_file, sep = "\t", header = TRUE, as.is = TRUE)
     } else { 
    # Throw an error if file cannot be found
    stop("Parameter table does not seem to exist")
    }

  # Identify the parameter name
  name_para <- colnames(param_table)

  # Name of parameter used as input in OpenMalaria: Remove scenario name, and add one parameter for each monthly EIR values
  param_names <- c(name_para[3:length(name_para)], paste0("seasonality", 1:12))

  # In case of the ACTS constrained analysis: Remove scenario name,  add one parameter for each monthly EIR values, and  add also a parameter for dosage
  if (pm$opts$Type_analysis == "Constrained" & pm$opts$drug == "ACT_SR") {
    param_names <- c(name_para[3:length(name_para)], "Dosage_S", paste0("seasonality", 1:12))
  }

  # Name of parameter used for the Global sensitivity analysis: Remove scenario name, but keep only one parameter for the seasonality profile
  param_names_2 <- c(name_para[3:length(name_para)], "seasonality")

  # Keep  the list of scenario
  
  # list of scenario without adjustment of the parameter values (table used for analysis later)
  Liste_scenario <- data.frame(matrix(ncol = length(param_names_2) + 1, nrow = nrow(param_table)))
  colnames(Liste_scenario) <- c(param_names_2, "scenario_name")

  # list of scenario with  adjustment of the parameter values for the input of OpenMalaria
  Liste_scenario_adjusted <- data.frame(matrix(ncol = length(param_names_2) + 1, nrow = nrow(param_table)))
  colnames(Liste_scenario_adjusted) <- c(param_names_2, "scenario_name")

  # Initiate progress bar
  pb <- txtProgressBar(min = 0, max = nrow(param_table), initial = 0, width = 100, style = 3)

  # Loop through table rows
  for (i in 1:nrow(param_table)) { # i<-1

    # select the ith scenario
    this_scen <- param_table[i, ]

    # select the seasonality profile
    this_season <- this_scen$seasonality

    #---- Creat the table used for the sensitivity analysis (without adjustment of parameter)----#
   
     # Insert the monthly EIR values
    seasonality_eir_2 <- pm$settings$seasonality[[this_season]] * this_scen$eir

    # Merge all parameter values without adjustment (for the output table)
    param_values_2 <- c(this_scen[3:length(name_para)])
    param_values_2 <- c(param_values_2, seasonality_eir_2)


    # ---- Creat the table of the parameter used for the input of OpenMalaria ---- #

    # adjust fitness values
    this_scen$Fitness <- adjust_Fitness(this_scen$Fitness) # see myfunction.R

    # adjust EIR
    if (this_season == "sesonality1") {
      this_scen$eir <- adjust_EIR(this_scen$eir, this_scen$Access) # see myfunction.R
    } else {
      this_scen$eir <- adjust_EIR_seasonality(this_scen$eir, this_scen$Access) # see myfunction.R
    }

    # Adjustement specific to long-acting drug
    if (pm$opts$drug == "long") {

      #  adjust IC50R
      this_scen$Resistance_Level <- this_scen$IC50_S * this_scen$Resistance_Level

      # adjuste half_life
      beta <- log(2) / this_scen$half_life
      a12 <- 8.46242774566474
      a21 <- 3.3058064516129035
      m <- 0.25
      W <- 60

      KK <- (beta^2 - a12 * beta - a21 * beta) / (beta - a21)
      this_scen$half_life <- KK * (W^m)


      # adjust the last treatment dosage (treatment adherence)
      if (pm$opts$Type_analysis == "Constrained") {
        this_scen$Dosage <- this_scen$Dosage * this_scen$Dosage_2
      }
    }

    # adjustment specific to long-acting drug
    if (pm$opts$drug == "short") {

      # adjust resistance level
      this_scen$Resistance_Level <- this_scen$MKR / this_scen$Resistance_Level

      # adjust half life
      this_scen$half_life <- log(2) / this_scen$half_life

      # adjust the last treatment dosage (treatment adherence)
      if (pm$opts$Type_analysis == "Constrained") {
        this_scen$Dosage <- this_scen$Dosage
      }
    }

    # adjustement specific to ACT
    if (pm$opts$drug == "ACT_SR") {

      # adjust resistance level of short acting drug
      this_scen$Resistance_Level <- this_scen$MKR_short / this_scen$Resistance_Level

      # adjust half life of short acting drug
      this_scen$half_life_short <- log(2) / this_scen$half_life_short

      # adjust the last treatment dosage of short acting drug(treatment adherence)
      if (pm$opts$Type_analysis == "Constrained") {
        this_scen$Dosage_S <- this_scen$Dosage * 4
      }

      # adjust resistance level of long acting drug
      this_scen$Resistance_Level_long <- this_scen$IC50_S_long * this_scen$Resistance_Level_long

      # adjust half_life of long acting drug
      beta <- log(2) / this_scen$half_life_long
      a12 <- 8.46242774566474
      a21 <- 3.3058064516129035
      m <- 0.25
      W <- 60
      KK <- (beta^2 - a12 * beta - a21 * beta) / (beta - a21)
      this_scen$half_life_long <- KK * (W^m)

      # adjust last treatment dosage of long acting drug (treatment adherence)
      if (pm$opts$Type_analysis == "Constrained") {
        this_scen$Dosage <- this_scen$Dosage * this_scen$Dosage_2
      }
    }

    # Multiple seasonality vector through by annual EIR value
    seasonality_eir <- pm$settings$seasonality[[this_season]] #* this_scen$eir

    # select all the parameter values use as input
    param_values <- c(this_scen[3:length(name_para)])
    param_values <- c(param_values, seasonality_eir)

    # select all the parameter values use as input when constrained and ACTs
    if (pm$opts$Type_analysis == "Constrained" & pm$opts$drug == "ACT_SR") {
      param_values <- c(this_scen[3:(length(name_para) + 1)]) # add +1 for the additional parameter for dosage
      param_values <- c(param_values, seasonality_eir)
    }

    # Format parameter values to not use scientific notation
    param_format <- format(param_values,
      scientific = FALSE,
      trim = TRUE, drop0trailing = TRUE)

    # Create sed command replacement pattern for each parameter
    sed_patterns <- paste0("s$@", param_names, "@$", param_format, "$g;")

    # save file
    file_name <- paste0(this_scen$scenario_name, "_", this_scen$seed, ".txt")
    file_path <- paste0(pm$pth$sim_sed, file_name)
    write.table(sed_patterns,
      file = file_path, quote = FALSE,
      col.names = FALSE, row.names = FALSE)

    # Conserve all parameter values of each scenario (with and without adjustment) in one table
    Liste_scenario[i, ] <- c(param_values_2[1:(length(param_values_2) - 12)], this_season, file_name)
    Liste_scenario_adjusted[i, ] <- c(param_values[1:(length(param_values) - 12)], this_season, file_name)

    # Update progress bar
    setTxtProgressBar(pb, i)
  } # Close seed loop

  # Save all info on the scenarios without adjustment of the parameter values
  file_name <- paste0(pm$pth$sim_sed, "scenarios_", ".txt")
  write.table(Liste_scenario, file = file_name, row.names = FALSE)

  # Save all info on the scenarios with adjustment of the parameter values
  file_name <- paste0(pm$pth$sim_sed, "scenarios_adjusted_", ".txt")
  write.table(Liste_scenario_adjusted, file = file_name, row.names = FALSE)

  # Close progress bar
  close(pb)

  return(nrow(param_table))
}

# -------------------------------------------------------------
# Generate scenario xml files from base file using sed patterns.
# -------------------------------------------------------------
generate_xml_files <- function(pm, n_jobs, file_path = NA,
                               sed_path = NA, xml_path = NA) {

  # message to the console
  message("  - Generating xml files")

  # Create a new log file for the cluster jobs
  log_file <- file.path(pm$pth$log_files, "scicore_log.txt")
  create_bash_log(log_file)

  # Add country sub directory to file paths
  sed_country <- paste0(sed_path)
  xml_country <- paste0(xml_path)

  # creat sequence of xml to be run at the same time
  n_seq <- 100
  n_submit <- ceil(n_jobs / n_seq)

  # message to console
  message("   ~ ", thou_sep(n_submit), " sets of ", thou_sep(n_seq))

  # Construct slurm array command for running in parallel
  slurm_array <- paste0("--array=1-", n_submit, "%", 300)

  # Concatenate system command
  sys_command <- paste(
    "sbatch", slurm_array, "bash_make_xmls.sh",
    n_seq, file_path, sed_country, xml_country, log_file)

  # Invoke this command
  system(sys_command)

  # Wait for all cluster jobs to complete
  wait_for_jobs(pm, log_file, n_jobs)
}
