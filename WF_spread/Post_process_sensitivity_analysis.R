#######################################################################################
# Post process results of the global sensitivity analysis                             #
#                                                                                     #
# Task : Reorganise the result of the sensitivity analysis from a list to a table     #
# Input: Output of the sensitivity analysis in a list format                          #
# output: Output of the sensitivity analysis in a table format                        #
#                                                                                     #
# Author:  thiery.masserey@swisstph.ch                                                #                                      ##
#######################################################################################

#---------------------------------------------------------------
# Function to organise the Sobol indices from a list to a table.
#---------------------------------------------------------------
Post_process_sensitivity <- function(Results_SA) {
  
  # Select the parameter name of the sensitivity analysis depending on drug modeled and if constrained or not
  if (pm$opts$Type_analysis == "Not_constrained") {
   pm$settings$Dosage <- pm$settings$Access <- pm$settings$eir <- pm$settings$Resistance_Level <- 1
   if (pm$opts$drug == "long") {
      parameters <- c("Fitness", "Access", "EIR", "Resistance level", "half-life", "MKR", "Diangostic", "Cmax/IC50")
    }
   if (pm$opts$drug == "ACT_SR") {
      parameters <- c("Fitness", "Access", "EIR", "Resistance level_long", "half-life_long", "MKR_long", "Diangostic", "IC50_S_short", "Resistance level_short", "half-life_short", "MKR_short", "Cmax/IC50")
    }
   if (pm$opts$drug == "short") {
      parameters <- pm$prog$Parameter
    }
  }
  if (pm$opts$Type_analysis == "Constrained") {
    if (pm$opts$drug == "long") {
      parameters <- c("Fitness", "half-life", "MKR", "Diangostic", "Cmax/IC50")
    }
    if (pm$opts$drug == "ACT_SR") {
      parameters <- c("Fitness", "half-life_long", "MKR_long", "Diangostic", "IC50_S_short", "half-life_short", "MKR_short", "Cmax/IC50")
    }
    if (pm$opts$drug == "short") {
      parameters <- pm$prog$Parameter
    }
  }

  # ---- Define the different variable of our ouput table and their length ----
  
  # name of the arm
  Setting_names <- rep(0, length(pm$settings$Access) * length(pm$settings$Resistance_Level) * length(parameters) * 2 * length(pm$settings$Dosage) * length(pm$settings$seasonality) * length(pm$settings$eir))

  # name of the parameter  
  factors <- rep(parameters, length(pm$settings$eir) * length(pm$settings$Access) * length(pm$settings$Resistance_Level) * 2 * length(pm$settings$Dosage) * length(pm$settings$seasonality))
  
  # if first/ total order indicies
  First <- rep(0, length(pm$settings$eir) * length(pm$settings$Access) * length(pm$settings$Resistance_Level) * length(parameters) * 2 * length(pm$settings$Dosage) * length(pm$settings$seasonality))
  
  # value of the sobol indices
  Effect <- rep(0, length(pm$settings$eir) * length(pm$settings$Access) * length(pm$settings$Resistance_Level) * length(parameters) * 2 * length(pm$settings$Dosage) * length(pm$settings$seasonality))
  
  # maximum values of the sobol indices
  MAX <- rep(0, length(pm$settings$eir) * length(pm$settings$Access) * length(pm$settings$Resistance_Level) * length(parameters) * 2 * length(pm$settings$Dosage) * length(pm$settings$seasonality))
  
  # minimum values of the sobol indices
  MIN <- rep(0, length(pm$settings$eir) * length(pm$settings$Access) * length(pm$settings$Resistance_Level) * length(parameters) * 2 * length(pm$settings$Dosage) * length(pm$settings$seasonality))
  
  # define the value for the factors that were constrained
  EIR <- rep(0, length(pm$settings$eir) * length(pm$settings$Access) * length(pm$settings$Resistance_Level) * length(parameters) * 2 * length(pm$settings$Dosage) * length(pm$settings$seasonality))
  Resistance_level <- rep(0, length(pm$settings$eir) * length(pm$settings$Access) * length(pm$settings$Resistance_Level) * length(parameters) * 2 * length(pm$settings$Dosage) * length(pm$settings$seasonality))
  Treatment_access <- rep(0, length(pm$settings$eir) * length(pm$settings$Access) * length(pm$settings$Resistance_Level) * length(parameters) * 2 * length(pm$settings$Dosage) * length(pm$settings$seasonality))
  Dosage <- rep(0, length(pm$settings$eir) * length(pm$settings$Access) * length(pm$settings$Resistance_Level) * length(parameters) * 2 * length(pm$settings$Dosage) * length(pm$settings$seasonality))
  Seasonality <- rep(0, length(pm$settings$eir) * length(pm$settings$Access) * length(pm$settings$Resistance_Level) * length(parameters) * 2 * length(pm$settings$Dosage) * length(pm$settings$seasonality))

  # creat the parameter table
  data <- cbind(Setting_names, as.character(factors), First, Effect, MAX, MIN, EIR, Resistance_level, Treatment_access, Dosage, Seasonality)
    
  # define the columns names
  colnames(data)[1:11] <- c("Setting", "Factor", "First", "Effect", "MAX", "MIN", "EIR", "Resistance_level", "Treatment_access", "Dosage", "Seasonality")
    
  #---- Restructure the table ----
    
  # define the number of arms (settings)pPost-Processed by the loop 
  n_settings <- 1
    
  # Do a loop across all arms
  for (this_season in names(pm$settings$seasonality)) {
    for (this_Dosage in pm$settings$Dosage) {
      for (this_Access in pm$settings$Access) {
        for (this_Resistance in pm$settings$Resistance_Level) { 
          for (this_eir in pm$settings$eir) {

            
            # Define the arm name
            if (pm$opts$Type_analysis == "Constrained") {
              setting_name <- paste(pm$opts$om_outcome, this_season, "EIR", this_eir, "Access", this_Access, "Resistance", this_Resistance, "Dosage", this_Dosage, sep = "_") ### remove eir paste(pm$opts$om_outcome, "eir", this_eir, this_season, sep = "_")
            } else {
              setting_name <- paste(pm$opts$om_outcome, this_season, sep = "_") ### remove eir paste(pm$opts$om_outcome, "eir", this_eir, this_season, sep = "_")
            }

            # Load the Sensitivity analysis results for this arms
            SA_file <- paste0(pm$pth$results, setting_name, "_gp.RData")
            SA <- readRDS(SA_file)
            
            # other options: 
            # SA<-Results_SA[setting_name]
            # SA<-list.ungroup(SA, level = 1L)
            
            # transfer the  data from the sensitivity analysis output into the data frame
            data[c(n_settings:(n_settings + length(parameters) * 2 - 1)), 1] <- setting_name
            data[n_settings:(n_settings + length(parameters) - 1), 3] <- SA$S$original[1:length(parameters)]
            data[(n_settings + length(parameters)):(n_settings + length(parameters) * 2 - 1), 3] <- SA$T$original[1:length(parameters)]
            data[c(n_settings:(n_settings + length(parameters) - 1)), 4] <- "First"
            data[c((n_settings + length(parameters)):(n_settings + length(parameters) * 2 - 1)), 4] <- "Total"
            data[n_settings:(n_settings + length(parameters) - 1), 5] <- SA$S$`max. c.i.`[1:length(parameters)]
            data[(n_settings + length(parameters)):(n_settings + length(parameters) * 2 - 1), 5] <- SA$T$`max. c.i.`[1:length(parameters)]
            data[n_settings:(n_settings + length(parameters) - 1), 6] <- SA$S$`min. c.i.`[1:length(parameters)]
            data[(n_settings + length(parameters)):(n_settings + length(parameters) * 2 - 1), 6] <- SA$T$`min. c.i.`[1:length(parameters)]
            data[c(n_settings:(n_settings + length(parameters) * 2 - 1)), 7] <- this_eir
            data[c(n_settings:(n_settings + length(parameters) * 2 - 1)), 8] <- this_Resistance
            data[c(n_settings:(n_settings + length(parameters) * 2 - 1)), 9] <- this_Access
            data[c(n_settings:(n_settings + length(parameters) * 2 - 1)), 10] <- this_Dosage
            data[c(n_settings:(n_settings + length(parameters) * 2 - 1)), 11] <- this_season
            
            # update the number of setting we went across with the post processing
            n_settings <- n_settings + length(parameters) * 2
          }
        }
      }
    }
  }

  # Define data into a data frame 
  data <- as.data.frame(data)
  
  # Define the column to be numeric
  data$First <- as.numeric(levels(data$First))[data$First]
  data$MAX <- as.numeric(levels(data$MAX))[data$MAX]
  data$MIN <- as.numeric(levels(data$MIN))[data$MIN]
  data$First[data$First <= 0] <- 0
  
  # save the Dataframe of Sobol indices for each arm and parameter
  name_file <- paste0(pm$pth$results, "Sobol_indices", ".txt")
  write.table(data,
    file = name_file, sep = "\t",
    quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  # return the  dataframe of sobol indices for each arm and parameter
  return(data)
  
}

#------------------------------------------- ---------------------------
# Function to organise the effect of each factor from a list to a table.
#-----------------------------------------------------------------------
Post_process_sensitivity_2 <- function(Results_SA) {
  
  # Select the parameter name of the sensitivity analysis depending on drug modeled and if constrained or not
  if (pm$opts$Type_analysis == "Not_constrained") {
    pm$settings$Dosage <- pm$settings$Access <- pm$settings$eir <- pm$settings$Resistance_Level <- 1
    
    if (pm$opts$drug == "long") {
      parameters <- c("Fitness", "Access", "EIR", "Resistance level", "half-life", "MKR", "Diangostic", "Cmax/IC50")
    }
    if (pm$opts$drug == "ACT_SR") {
      parameters <- c("Fitness", "Access", "EIR", "Resistance level_long", "half-life_long", "MKR_long", "Diangostic", "IC50_S_short", "Resistance level_short", "half-life_short", "MKR_short", "Cmax/IC50")
    }
    if (pm$opts$drug == "short") {
      parameters <- pm$prog$Parameter
    }
  }
  if (pm$opts$Type_analysis == "Constrained") {
    if (pm$opts$drug == "long") {
      parameters <- c("Fitness", "half-life", "MKR", "Diangostic", "Cmax/IC50")
    }
    if (pm$opts$drug == "short") {
      parameters <- pm$prog$Parameter
    }
    if (pm$opts$drug == "ACT_SR") {
      parameters <- c("Fitness", "half-life_long", "MKR_long", "Diangostic", "IC50_S_short", "half-life_short", "MKR_short", "Cmax/IC50")
    }
  }

  
  # Creat a parameter table that will contains the estimate selection coefficient over the parameter range for each arm
  Quantil_final_final <- matrix(NA, nrow = 0, ncol = 5)
  
  # give name to the different columns
  colnames(Quantil_final_final) <- c("L", "M", "U", "x", "G")

  # Loop across each arm
  for (this_season in names(pm$settings$seasonality)) {
    for (this_Dosage in pm$settings$Dosage) {
      for (this_Access in pm$settings$Access) {
        for (this_Resistance in pm$settings$Resistance_Level) { 
          for (this_eir in pm$settings$eir) { 

            # define the name of the arm
            if (pm$opts$Type_analysis == "Constrained") {
              setting_name <- paste(pm$opts$om_outcome, this_season, "EIR", this_eir, "Access", this_Access, "Resistance", this_Resistance, "Dosage", this_Dosage, sep = "_") ### remove eir paste(pm$opts$om_outcome, "eir", this_eir, this_season, sep = "_")
            } else {
              setting_name <- paste(pm$opts$om_outcome, this_season, sep = "_") ### remove eir paste(pm$opts$om_outcome, "eir", this_eir, this_season, sep = "_")
            }
            
            # load the results of the sensitivity analyis for this arm
            SA_file <- paste0(pm$pth$results, setting_name, "_gp.RData")
            SA <- readRDS(SA_file)
            
            # other options:
            # SA<-Results_SA[setting_name]
            # SA<-list.ungroup(SA, level = 1L)

            # Select the predicted selection coefficient into a vector
            Y <- SA$y
            
            # select the input parameter into a table
            X <- matrix(NA, ncol = length(SA$X), nrow = length(Y))
            colnames(X) <- colnames(SA$X)
            for (i in 1:length(SA$X)) { 
              X[, i] <- SA$X[, i]
            }

            # For each parameter, divide the parameter into N fraction, and esitmate the mean selection coefficient in this fraction
            
            # Define the number of fraction
            N <- 40
            
            # Creat a vector that will contain the boundary of fraciton of the parameter space + median selection coeffcient 
            Quantil_final <- matrix(NA, nrow = 0, ncol = 5)
            
            
            # loop across each parameters
            for (i in 1:length(SA$X)) {
              
              # Creat a datafram of parameter value and  selection coefficient
              DF <- as.data.frame(cbind(X[, i], Y))
              
              # cut the parameter space into N piece
              splitFitness <- seq(min(X[, i]), max(X[, i]), length = N)
              
              # Creat a table for this specific parameter
              Quantiles_1 <- matrix(NA, nrow = N, ncol = 3)

              # estimate the median selection coeffcient in each fraction of the parameter space
              for (k1 in 1:(N - 1)) {
                thisXY <- filter(DF, X[, i] > splitFitness[k1] & X[, i] < splitFitness[k1 + 1])
                Quantiles_1[k1, ] <- quantile(thisXY$Y, c(0.25, 0.5, 0.75))
              }

              # Transform data into dataframe
              Quantiles_1 <- as.data.frame(cbind(Quantiles_1, 1:40))
              
              # Select the parameter name
              Quantiles_1$G <- colnames(X)[i]
              
              # Add columns name
              colnames(Quantiles_1) <- c("L", "M", "U", "x", "G")

              # Update the table to contain the data of all parameters
              Quantil_final <- rbind(Quantil_final, Quantiles_1)
            
            }
            
            # Add information about the arm
            Quantil_final$EIR <- this_eir
            Quantil_final$Resistance_level <- this_Resistance
            Quantil_final$Access <- this_Access
            Quantil_final$Dosage <- this_Dosage
            Quantil_final$Seasonality <- this_season
            
            # Update the table to contain the data of all arm
            Quantil_final_final <- rbind(Quantil_final_final, Quantil_final)
          }
        }
      }
    }
  }

  # Save the table of mean seleciton coeffcient over the parameter range of each paramter, for each arm
  name_file <- paste0(pm$pth$results, "Factors_effect", ".txt")
  write.table(Quantil_final_final,
    file = name_file, sep = "\t",
    quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  # Return the table of mean selection coeffcient over the parameter range of each paramter, for each arm
  return(Quantil_final_final)
}