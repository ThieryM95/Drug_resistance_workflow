###################################################################################
# Perform the sensitivity analysis                                                #
#                                                                                 #
# Task : Perfom the global sensitivity analysis of the GP                         #
# Input: fitted GP                                                                #
# Output: sobol indices + direction of effect                                     #
#                                                                                 #
# authors: thiery.masserey@swisstph.ch adapted from monica.golumbeanu@swisstph.ch #                                      #
###################################################################################

# load the package
library("sensitivity")
library("multisensi")
library("lhs")
library("hetGP")
library("tgp")
library("dplyr")
library("ggplot2")
library("reshape2")
library("gridExtra")
library("hrbrthemes")
library("grid")
library("gridExtra")
library("rlist")

#--------------------------------------------------------------------
# function to preapre data for sensitivity analysis and save results.
#--------------------------------------------------------------------
sensitivity_analysis <- function(pm) {

  # creat a liste to store the results of the SA
  indices <- list()

  # define the parameter range
  param_ranges <- cbind(pm$prog$min, pm$prog$max)
  if (pm$opts$Type_analysis == "Constrained" & pm$opts$drug == "ACT_SR") {
    param_ranges[8, 1] <- 0.08
  }

  # define the parameter names
  row.names(param_ranges) <- pm$prog$prog_names

  # define the column names
  colnames(param_ranges) <- c("min", "max")

  # number of sampled points
  sa_n <- pm$opts$sa_n

  # loop across each arms
  if (pm$opts$Type_analysis == "Not_constrained") {
    pm$settings$Dosage <- pm$settings$Access <- pm$settings$eir <- pm$settings$Resistance_Level <- pm$settings$Resistance_Level_long <- 1
  }
  if (pm$opts$drug == "short" | pm$opts$drug == "long") {
    pm$settings$Resistance_Level_long <- 1
  }

  for (this_Dosage in pm$settings$Dosage) {
    for (this_Access in 0.5) { 
      for (this_Resistance in pm$settings$Resistance_Level) { 
        for (this_eir in pm$settings$eir) { 
          for (this_season in names(pm$settings$seasonality)) { 
            for (this_Resistance_L in (pm$settings$Resistance_Level_long)) {

              # Define the name of the arm
              if (pm$opts$Type_analysis == "Constrained") {
                setting_name <- paste(pm$opts$om_outcome, this_season, "EIR", this_eir, "Access", this_Access, "Resistance", this_Resistance, "Dosage", this_Dosage, sep = "_") ### remove eir paste(pm$opts$om_outcome, "eir", this_eir, this_season, sep = "_")
              } else {
                setting_name <- paste(pm$opts$om_outcome, this_season, sep = "_") ### remove eir paste(pm$opts$om_outcome, "eir", this_eir, this_season, sep = "_")
              }

              # Message to the console
              message("  - ", setting_name)

              # Load the GP of the define arm
              gp_file <- paste0(pm$pth$gp_models, setting_name, "_gp.RData")
              trained_model <- readRDS(gp_file)
              trained_model <- trained_model$gp_model

              # Do the sensitivity analysis
              sobol_indices <- calc_sobol_idx(trained_model, param_ranges, sa_n, this_Resistance, this_Dosage, pm$opts$Type_analysis, pm$opts$drug) # see function bellow

              # Save the results in a dataframe
              SA_file <- paste0(pm$pth$results, setting_name, "_gp.RData")
              saveRDS(sobol_indices, file = SA_file)

              # Save the results as an output of this function
              indices$setting_name <- sobol_indices
              names(indices)[names(indices) == "setting_name"] <- setting_name
            }
          }
        }
      }
    }
  }
  
  # return the sobole indice
  return(indices)
}

#----------------------------------------------
# function to perform the sensitivity analysis.
#----------------------------------------------
calc_sobol_idx <- function(trained_model, param_ranges, sa_n, this_Resistance, this_Dosage, Type_analysis, drug) {

  # wrapper function for the GP_model prediction
  GP_f <- function(X) {
    out <- predict(x = as.matrix(X), object = trained_model)
    return(out$mean)
  }

  # construct the two random lhs samples
  X1 <- as.data.frame(lhs(sa_n, as.matrix(param_ranges)))
  X2 <- as.data.frame(lhs(sa_n, as.matrix(param_ranges)))

  # define the columns name
  colnames(X1) <- pm$prog$Parameter
  colnames(X2) <- pm$prog$Parameter

  # if long acting drug  or ACT estimate Cmax_EC50 for the long actiong drug based on dosage and treatment adherence (treatment adherence only if constrain analysis)
  if (Type_analysis == "Constrained" & drug == "long") {
    
    # estimate Cmax_EC50 for each dataset if a 60% treatment adhenrence
    if (this_Dosage == 0) {
      X1$C_max_IC50 <- ((-8.368621e-17) + X1$Dosage_2 * (4.048628e-03)) / (X1$IC50_S)
      X2$C_max_IC50 <- ((-8.368621e-17) + X2$Dosage_2 * (4.048628e-03)) / (X2$IC50_S)
    
    # estimate Cmax_EC50 for each dataset if a 100% treatment adhenrence
    } else {
      X1$C_max_IC50 <- ((-1.004235e-16) + X1$Dosage_2 * (5.440722e-03)) / (X1$IC50_S)
      X2$C_max_IC50 <- ((-1.004235e-16) + X2$Dosage_2 * (5.440722e-03)) / (X2$IC50_S)
    }
    
    # Select the columns of interest (Cmax insteand of Ec50 and dosage)
    X1 <- X1[, c("Fitness", "half_life", "MKR", "Diagnostic", "C_max_IC50")]
    X2 <- X2[, c("Fitness", "half_life", "MKR", "Diagnostic", "C_max_IC50")]
  }

  if (Type_analysis == "Not_constrained" & drug == "long") {

    # estimate Cmax_EC50 for each dataset with a 100% treatment adhenrence
    X1$C_max_IC50 <- ((-1.004235e-16) + X1$Dosage_2 * (5.440722e-03)) / (X1$IC50_S)
    X2$C_max_IC50 <- ((-1.004235e-16) + X2$Dosage_2 * (5.440722e-03)) / (X2$IC50_S)

    # Select the columns of interest (Cmax insteand of Ec50 and dosage)
    X1 <- X1[, c("Fitness", "Access", "eir", "Resistance_Level", "half_life", "MKR", "Diagnostic", "C_max_IC50")]
    X2 <- X2[, c("Fitness", "Access", "eir", "Resistance_Level", "half_life", "MKR", "Diagnostic", "C_max_IC50")]
  }

  if (Type_analysis == "Not_constrained" & drug == "ACT_SR") {
    
    # estimate Cmax_EC50 for each dataset with a 100% treatment adhenrence
    X1$C_max_IC50 <- ((-1.004235e-16) + X1$Dosage_2_long * (5.440722e-03)) / (X1$IC50_S_long)
    X2$C_max_IC50 <- ((-1.004235e-16) + X2$Dosage_2_long * (5.440722e-03)) / (X2$IC50_S_long)

    # Select the columns of interest (Cmax insteand of Ec50 and dosage)
    X1 <- X1[, c("Fitness", "Access", "eir", "Resistance_Level_long", "half_life_long", "MKR_long", "Diagnostic", "IC50_S_short", "Resistance_Level", "half_life_short", "MKR_short", "C_max_IC50")]
    X2 <- X2[, c("Fitness", "Access", "eir", "Resistance_Level_long", "half_life_long", "MKR_long", "Diagnostic", "IC50_S_short", "Resistance_Level", "half_life_short", "MKR_short", "C_max_IC50")]
  }

  if (Type_analysis == "Constrained" & drug == "ACT_SR") {
    
    # estimate Cmax_EC50 for each dataset if a 60% treatment adhenrence
    if (this_Dosage == 0) {
      X1$C_max_IC50 <- ((-8.368621e-17) + X1$Dosage_2_long * (4.048628e-03)) / (X1$IC50_S_long * pm$settings$Resistance_Level_long)
      X2$C_max_IC50 <- ((-8.368621e-17) + X2$Dosage_2_long * (4.048628e-03)) / (X2$IC50_S_long * pm$settings$Resistance_Level_long)
    
    # estimate Cmax_EC50 for each dataset if a 100% treatment adhenrence
    } else {
      X1$C_max_IC50 <- ((-1.004235e-16) + X1$Dosage_2_long * (5.440722e-03)) / (X1$IC50_S_long * pm$settings$Resistance_Level_long)
      X2$C_max_IC50 <- ((-1.004235e-16) + X2$Dosage_2_long * (5.440722e-03)) / (X2$IC50_S_long * pm$settings$Resistance_Level_long)
    }

    # Select the columns of interest (Cmax insteand of Ec50 and dosage)
    X1 <- X1[, c("Fitness", "half_life_long", "MKR_long", "Diagnostic", "IC50_S_short", "half_life_short", "MKR_short", "C_max_IC50")]
    X2 <- X2[, c("Fitness", "half_life_long", "MKR_long", "Diagnostic", "IC50_S_short", "half_life_short", "MKR_short", "C_max_IC50")]
  }

  # estimate the Sobol indices and direction of effect using sobol jasen
  SA <- soboljansen(model = GP_f, as.data.frame(X1), as.data.frame(X2), nboot = 15000)
  
  # return sobol sensitivity indices + direction of effect
  return(SA)
}
