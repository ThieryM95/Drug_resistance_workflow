#######################################################################
# Post process of the GP prediction for a better visulaisaiton        #
#                                                                     #
# Task : Reorganise the gp accuracy result from a list to a dataframe #
# Input: GP prediction and true value from  the train dataset         #
# Output: Table of the GP test and prediction for each arm            #
#                                                                     #
# Author:  thiery.masserey@swisstph.ch                                #
#######################################################################

# load package
library(SummarizedExperiment)

#--------------------------------------------------------------------------
# function to Reorganise the gp accuracy result from a list to a dataframe.
#--------------------------------------------------------------------------
Post_process_GP <- function(Results_gp) {

  # Creat the dataframe
  precision_final <- precision_3 <- precision_2 <- matrix(NA, nrow = 0, ncol = 8)

  # Add the columns name
  colnames(precision_2)[1:8] <- c("Test_True", "Test_predicted", "iteration", "eir", "seasonality", "dosage", "access", "resistance_level")

  # Define some parameter too one for non constrained analysis
  if (pm$opts$Type_analysis == "Not_constrained") {
    pm$settings$Dosage <- pm$settings$Access <- pm$settings$eir <- pm$settings$Resistance_Level <- 1
  }
  
  # Loop across the different  arms
  for (iter in 1:length(Results_gp)) {
    
    # select data of the ith iteration
    z <- Results_gp[[iter]]
    
    # Loop across the different  arms
    for (this_season in names(pm$settings$seasonality)) {
      for (this_Dosage in pm$settings$Dosage) {
        for (this_Access in pm$settings$Access) {
          for (this_Resistance in pm$settings$Resistance_Level) {
            for (this_eir in pm$settings$eir) {

              # define the arm name
              if (pm$opts$Type_analysis == "Constrained") {
                setting_name <- paste(pm$opts$om_outcome, this_season, "EIR", this_eir, "Access", this_Access, "Resistance", this_Resistance, "Dosage", this_Dosage, sep = "_") ### remove eir paste(pm$opts$om_outcome, "eir", this_eir, this_season, sep = "_")
              } else {
                setting_name <- paste(pm$opts$om_outcome, this_season, sep = "_") ### remove eir paste(pm$opts$om_outcome, "eir", this_eir, this_season, sep = "_")
              }

              # select the data from the arm
              precision <- z[setting_name]

              # Unlist the results
              precision <- list.ungroup(precision, level = 1L)

              # select the variable true value and predicted value
              Test_True <- precision$actual_test
              Test_predicted <- precision$predict_test
              iteration <- iter
              eir <- this_eir
              seasonality <- this_season
              dosage <- this_Dosage
              access <- this_Access
              resistance_level <- this_Resistance

              # creat the data  with all variable
              precision_2 <- data.frame(Test_True, Test_predicted, iteration, eir, seasonality, dosage, access, resistance_level)

              # creat a data that contain all the data from all the arms
              precision_3 <- rbind(precision_3, precision_2)
            }
          }
        }
      }
    }
  }

  # save the datas
  name_file <- paste0(pm$pth$gp_models, "Precision", ".txt")
  write.table(precision_3,
    file = name_file, sep = "\t",
    quote = FALSE, col.names = TRUE, row.names = FALSE)

  # return the dataset
  return(precision_final)
}
