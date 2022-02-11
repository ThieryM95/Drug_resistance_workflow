######################################################################################
#      Generate the parameter tables for the simulation that will be run             #
#                                                                                    #
#                                                                                    #
# Task: This function create the table of input parameter for each simulation        #
#                                                                                    #  
# Input: Number of simulations, number of seeds, constrained parameter values,       #
#        parameter range                                                             #
#                                                                                    #
# Output: A table of the parameter value of each parameter for each simulation       #
#                                                                                    #
# authors: thiery.masserey@swisstph.ch adapted from andrewjames.shattock@swisstph.ch #                                      #
######################################################################################

# --------------------------------------------
# Generate the input table for the simulation.
# --------------------------------------------
generate_param_table <- function(pm, param_table = NULL) {

  # Message in console
  message("* Generating parameter table")

  # File name for this parameter table (see adaptive_sampling.R)
  param_table_name <- param_set_name(sample_num = pm$sample_num)
  param_file <- paste0(pm$pth$param_table, param_table_name, ".txt")

  # If sample number is == 0 the parameter table is generate via Latin Hypercube sampling
  if (pm$sample_num == 0) {

    # Continue if file doesn't exist OR we want to overwrite
    if (!file.exists(param_file) | pm$opts$do_resample == TRUE) {

      # ---- Generate latin hypercube samples for intervention coverage ----

      # Load the define parameter  ranges into a matrix
      cov_ranges <- matrix(c(pm$prog$min, pm$prog$max), ncol = 2)

      # Generate latin hypercube sample between the coverage limits
      param_table <- as.data.frame(lhs(pm$opts$lhc_samples, cov_ranges))

      # Variable names of parameters must match @parameter@ fields in base xml
      colnames(param_table) <- paste(pm$prog$prog_names)

      # ---- Define the parameter that were constrained (the different arms) ----

      # If not constrained analysis
      if (pm$opts$Type_analysis == "Not_constrained") {

        # Define the level and value for each parameter
        season <- pm$settings$seasonality
        
        # Merge into full factorial table
        setting_table <- data.frame(names(season))
        colnames(setting_table) <- c("seasonality")
      }

      # If constrained analysis and short or long acting drug
      if (pm$opts$Type_analysis == "Constrained" & pm$opts$drug != "ACT_SR") {
        
        # Define the level and value for each parameter
        eir <- pm$settings$eir
        season <- pm$settings$seasonality
        Resistance_Level <- pm$settings$Resistance_Level
        Access <- pm$settings$Access
        Dosage <- pm$settings$Dosage

        # Merge into full factorial table
        setting_table <- Reduce(merge, list(names(season), as.data.frame(eir), as.data.frame(Access), as.data.frame(Resistance_Level), as.data.frame(Dosage)))
        colnames(setting_table) <- c("seasonality", "eir", "Access", "Resistance_Level", "Dosage") # remove eir before c("eir", "seasonality") ######################################
      }

      # if constrained analysis and ACT
      if (pm$opts$Type_analysis == "Constrained" & pm$opts$drug == "ACT_SR") {
        
        # Define the level and value for each parameter
        eir <- pm$settings$eir
        season <- pm$settings$seasonality
        Resistance_Level_long <- pm$settings$Resistance_Level_long
        Resistance_Level <- pm$settings$Resistance_Level
        Access <- pm$settings$Access
        Dosage <- pm$settings$Dosage

        # Merge into full factorial table
        setting_table <- Reduce(merge, list(names(season), as.data.frame(eir), as.data.frame(Access), as.data.frame(Resistance_Level_long), as.data.frame(Resistance_Level), as.data.frame(Dosage)))
        colnames(setting_table) <- c("seasonality", "eir", "Access", "Resistance_Level_long", "Resistance_Level", "Dosage") # remove eir before c("eir", "seasonality") ######################################
      }

      # Merge the setting (arms) table and parameter table
      param_table <- merge(setting_table, param_table)

      # Incorporate a scenario name column
      scenario_name <- paste("scenario", 1:nrow(param_table), sep = "_")
      param_table <- cbind(scenario_name, param_table)
    }

    # Repeat the scenarios for a predefined number of seeds
    seed <- 1:pm$opts$n_seeds
    param_table <- merge(param_table, as.data.frame(seed))

    # Write table to file
    write.table(param_table,
      file = param_file, sep = "\t",
      quote = FALSE, col.names = TRUE, row.names = FALSE)

  } else {
    # Case where sample number is greater than 0 the parameter table is generate via adaptive sampling
    
    # Write extended parameter table to a new file
    write.table(param_table,
      file = param_file, sep = "\t",
      quote = FALSE, col.names = TRUE, row.names = FALSE
    )
  }
}
