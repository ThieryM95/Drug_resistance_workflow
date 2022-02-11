##################################################################
#     Estimate the selection coeffcient for each simulation      #
#                                                                #
# Task: Calculate the selection coefficient for each simulation  #
# Input: Post processed data output of OpenMalaria               #
# Output: Table of parameter input of each simulation and        #
#         the estimated selection coefficient                    #
#                                                                #
# Important: need to define here the number of time step         #
# of the survey that was used                                    #
# authors: thiery.masserey@swisstph.ch                           #                                      
##################################################################

# load function
library(grDevices)
library(pkr)

# ------------------------------------------------------
# Estimate the selection coefficient for each simulation.
# ------------------------------------------------------
SummaryResults <- function(pm) {

  # Message in the console
  message("  - Calulate the spread or esthablihsment ")

  # Load the list of scenairos
  Scenario_liste <- read.table(file.path(pm$pth$sim_sed, "scenarios_.txt"), header = T)

  # Load the list of outputs
  Output <- list.files(pm$pth$processed)

  # Define the time step
  time_step <- 30

  # Define the number of survey per parasite generation
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days

  # Define the number of survey per years
  Number_survey_years <- 365 / time_step # how many survey per years

  # Creat a new variable that will store the estimated value for the selection coefficient
  Scenario_liste$Indicator <- NA

  # Initiate progress bar
  pb <- txtProgressBar(
    min = 0, max = length(Output),
    initial = 0, width = 100, style = 3
  )

  # Start the loop across all the scenario
  for (i in 1:length(Output)) {

    # Download the data
    Output_data_name <- Output[i]
    Output_data_file_path <- file.path(pm$pth$processed, Output_data_name)
    Output_data <- read.table(Output_data_file_path, sep = ";")
    name_scenario <- gsub("PostProcess_", "", Output_data_name)
    name_scenario <- gsub("_out", "", name_scenario)
    
    # Estimate the rate of spread in setting without seasonality
    if (Scenario_liste$seasonality[Scenario_liste$scenario_name == name_scenario] == "sesonality1") {
      
      Scenario_liste$Indicator[Scenario_liste$scenario_name == name_scenario] <- spread_R2(Output_data) # see function bellow
    
    # Estimate the rate of spread in setting with seasonality (use a moving average)
    } else {
      Scenario_liste$Indicator[Scenario_liste$scenario_name == name_scenario] <- spread_R5(Output_data) # see function bellow
    }

    # update the progress bar
    setTxtProgressBar(pb, i)
    
    } # end the loop

  # Close the progress bar
  close(pb)

  # Save the  dataset in the summarized folder
  name_sumarised <- param_set_name(sample_num = pm$sample_num)
  summary_file <- paste0(name_sumarised, "_", ".txt")
  summary_path <- paste0(pm$pth$summarised, summary_file)
  write.table(Scenario_liste,
    file = summary_path, sep = "\t",
    col.names = TRUE, row.names = FALSE, quote = FALSE)

  # Merge the summarized results of each iterations
  full_table_name <- param_set_name(all_samples = TRUE)

  # Construct file name and path
  full_summary_file <- paste0(full_table_name, "_", ".txt")
  full_summary_path <- paste0(pm$pth$summarised, full_summary_file)

  # Initiate 'all sample' table
  full_param_table <- NULL

  # Process only necessary for adaptive sampling summary
  if (pm$sample_num > 0) {

    # Loop through previously generated samples
    for (i in (1:pm$sample_num - 1)) {

      # Construct sample results summary file name and path
      sample_file <- paste0(
        param_set_name(sample_num = i),
        "_", ".txt")

      # Load up the parameter table associated with this sample number
      sample_param_table <- read.table(paste0(pm$pth$summarised, sample_file),
        sep = "\t", header = TRUE, as.is = TRUE)
      
      # Colnames(sample_param_table)[11]<-"Resistance_Level_short"
      full_param_table <- rbind(full_param_table, sample_param_table)
    }
  }

  # Concatenate the param table we have just summarized
  full_param_table <- rbind(full_param_table, Scenario_liste)

  # Overwrite scenario names in this 'all samples' file
  full_param_table$scenario_name <- paste0("scenario_", 1:nrow(full_param_table))

  # (Over)-write full result tables across all samples
  write.table(full_param_table, full_summary_path,
    sep = "\t",
    col.names = TRUE, row.names = FALSE, quote = FALSE)

} # close function


# ------------------------------------------------------------------
# Function to estimate the selection coefficient in perinal setting.
# ------------------------------------------------------------------
spread_R2 <- function(Output_data) {

  # Define time step
  time_step <- 30

  # Number of survey per parasite generation
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days

  # Number of survey per years
  Number_survey_years <- 365 / time_step

  # Time step at which drug resistance is introduced
  Time_HS_Change <- round(30 * Number_survey_years)

  # Define time at which regression will start (when drug concentration of the previous drug is equal to zero + 1 parasite generations)
  time_start <- Output_data$Survey[Output_data$nHostDrugConcNonZero_2 == 0]
  time_start <- time_start[time_start >= (Time_HS_Change + Number_survey_generation)][1]

  # Define time at which regression end (2 years after the beginning of the regression)
  time_end <- time_start + 2 * Number_survey_years + 2

  # Select the data within this boundary (time step and frequency of resistant parasite)
  time_spread <- Output_data$Survey[Output_data$Survey > time_start & Output_data$Survey < time_end] / Number_survey_years # time spread is converted in years
  Measurment_R <- Output_data$Inoculation_R[Output_data$Survey > time_start & Output_data$Survey < time_end]

  # Remove data for which the frequency of the resistant genotype is to high or to low
  Measurment_R <- Measurment_R[Measurment_R <= 0.9 & Measurment_R >= 0.3]
  time_spread <- time_spread[1:length(Measurment_R)]

  # If there is at least two measurements to calculate the selection coefficient
  if (length(Measurment_R) >= 2 & is.na(sum(Measurment_R)) == F) {

    # Perform the selection coefficient
    slope <- lm(log(Measurment_R / (1 - Measurment_R)) ~ time_spread)$coefficients[2]
    slope <- slope / 6
    } else {

    # If there is not at least two measurements return NA
    slope <- NA
    }
  
  # return selection coefficient
  return(slope)
  
} # close function

# ------------------------------------------------------------------
# Function to estimate the selection coefficient in seasonal setting.
# ------------------------------------------------------------------
spread_R5 <- function(Output_data) { 

  # Define the number of time step
  time_step <- 30

  # Define the number of survey per generation
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days

  # Define the number of survey per years
  Number_survey_years <- 365 / time_step

  # define the time at which drug resistance is introduce
  Time_HS_Change <- round(30 * Number_survey_years)

  # Define time at which regression will start (when drug concentration of the previous drug is equal to zero + 1 parasite generations)
  time_start <- Output_data$Survey[Output_data$nHostDrugConcNonZero_2 == 0]
  time_start <- time_start[time_start >= (Time_HS_Change + Number_survey_generation)][1]

  # Define time at which regression will stop (2 years after the beginning of the regression + 6 more for the moving average)
  time_end <- time_start + 2 * Number_survey_years + 6

  # select the data within this boundary
  time_spread <- Output_data$Survey[Output_data$Survey > time_start & Output_data$Survey < time_end] / Number_survey_years # time spread is converted in years
  Measurment_R <- Output_data$Inoculation_R[Output_data$Survey > time_start & Output_data$Survey < time_end]
  
  # estimate the moving average (over 1 year period) of the frequency of the resistant genotype
  Measurment_MA <- 0
  time_spread_MA <- 0
  for (i in (1 + 6):(24)) {
    Measurment_MA[i - 6] <- mean(Measurment_R[(i - 6):(i + 6)])
    time_spread_MA[i - 6] <- time_spread[i]
  }
  
  plot(time_spread-30, Measurment_R, pch =19, xlab="Time since the end of the burn-in phase (years)", ylab="The logit of the relative frequency of\n the resistant genotype in inoculations", cex.lab=1.5,cex.axis=1.4, lwd=3,font.lab=2)
  points(time_spread_MA-30, Measurment_MA, col="blue", pch=19, lwd=3)
  box(lwd = 3)
  

  # check that the moving average have not to low frequency or to high frequency but include that need at maximum of 12 measurement to not bias the estimate
  Measurment_MA <- Measurment_MA[1:12]
  time_spread_MA <- time_spread_MA[1:12]
  Measurment_MA <- Measurment_MA[Measurment_MA >= 0.3 & Measurment_MA <= 0.99] #
  time_spread_MA <- time_spread_MA[1:length(Measurment_MA)]

  # if there is at least two measurements calculate the selection coefficient
  if (length(Measurment_MA) >= 2 & is.na(sum(Measurment_MA)) == F) {
    # perform the logistic regression
    slope <- lm(log(Measurment_MA / (1 - Measurment_MA)) ~ time_spread_MA)$coefficients[2]
    slope <- slope / 6
  } else {
    # if there is not at least two measurements return NA
    slope <- NA
  }
  
  # return selection coefficient
  return(slope)
  
} # close the loop

