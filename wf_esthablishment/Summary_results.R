#####################################################################
# Estimate the  probability of establishment for each simulation    #
#                                                                   #
# Task: Calculate the probability of establishment across seed      #
#                                                                   #
# Input: Post-processed data output of OpenMalaria                  #
#                                                                   #
# Output: Table of parameter input of each simulation and           #
#         estimated the probability of establishment                #
#                                                                   #
# author: thiery.masserey@swisstph.ch                               #
#####################################################################

# load function
library(grDevices)
library(pkr)

# ------------------------------------------------------
# Function to estimate the probability of establishment.
# ------------------------------------------------------
SummaryResults <- function(pm) {
  
  # Message in the console
  message("  - Calulate the esthablihsment ")

  #---- Estimate how many mutation get extinct in each simulation ----#

  # Load the input parameter table with list of scenarion
  Scenario_liste <- read.table(file.path(pm$pth$sim_sed, "scenarios_.txt"), header = T)

  # Load the list of output simulation
  Output <- list.files(pm$pth$processed)

  # Define the time setp of the survey
  time_step <- 30
  Number_survey_generation <- 60 / time_step
  Number_survey_years <- 365 / time_step

  # List of new variable that will store the info needed for each simulation
  Scenario_liste$eir_sim <- 0
  Scenario_liste$prevalance <- 0
  Scenario_liste$esthablisment <- NA
  Scenario_liste$Number_mutation <- NA
  Scenario_liste$Number_mutation_2 <- NA
  Scenario_liste$Indicator_2 <- 0

  # Initiate progress bar
  pb <- txtProgressBar(
    min = 0, max = length(Output),
    initial = 0, width = 100, style = 3
  )

  # Start a loop across each simulation
  for (i in 1:length(Output)) {

    # Load the data
    Output_data_name <- Output[i]
    Output_data_file_path <- file.path(pm$pth$processed, Output_data_name)
    Output_data <- read.table(Output_data_file_path, sep = ";")
    name_scenario <- gsub("PostProcess_", "", Output_data_name)
    name_scenario <- gsub("_out", "", name_scenario)

    # Estimate the number of mutation that emerge in the simulation using the formula from Hasting et al. (2020)
    Scenario_liste$Number_mutation[Scenario_liste$scenario_name == name_scenario] <- Number_mutation(Output_data) # see function bellow

    # Estimate the number of mutation that emerge in the simulation by compting the number of time we see a mutation occure
    Scenario_liste$Number_mutation_2[Scenario_liste$scenario_name == name_scenario] <- Number_mutation_2(Output_data) # see function bellow

    # Estimate if one mutation established
    Scenario_liste$esthablisment[Scenario_liste$scenario_name == name_scenario] <- Esthablishment(Output_data) # see function bellow

    # Estimate the rate of spread of the mutation that established
    Scenario_liste$Indicator_2[Scenario_liste$scenario_name == name_scenario] <- spread(Output_data) # see function bellow

    # Estimate simulated EIR and the prevalence
    Scenario_liste$eir_sim[Scenario_liste$scenario_name == name_scenario] <- mean(Output_data$simulatedEIR[(50 * Number_survey_years):(200 * Number_survey_years)]) * 73
    Scenario_liste$prevalance[Scenario_liste$scenario_name == name_scenario] <- mean(Output_data$nInfect[(50 * Number_survey_years):(200 * Number_survey_years)]) / 10000

    # Update the progress bar
    setTxtProgressBar(pb, i)
    
  }

  # close the progress bar
  close(pb)

  # Estimate the mean selection coefficient across seed to see if it was as expected
  for (this_s in Scenario_liste$Indicator) {
    Scenario_liste$mean_selection[Scenario_liste$Indicator == this_s] <- mean(Scenario_liste$Indicator_2[Scenario_liste$Indicator == this_s], na.rm = T)
  }

  # estimate the mean prevalence across seed to see if it was as expected
  for (this_s in Scenario_liste$Indicator) {
    Scenario_liste$mean_prevalance[Scenario_liste$Indicator == this_s] <- mean(Scenario_liste$prevalance[Scenario_liste$Indicator == this_s], na.rm = T)
  }

  # Save the  dataset in the sumarised folder
  name_sumarised <- param_set_name(sample_num = 0)
  summary_file <- paste0(name_sumarised, "_new", ".txt")
  summary_path <- paste0(pm$pth$summarised, summary_file)
  write.table(Scenario_liste,
    file = summary_path, sep = "\t",
    col.names = TRUE, row.names = FALSE, quote = FALSE
  )

  #---- Estimate the probability of establishment across seed ----#

  # Select only the columns of interest
  Scenario_liste_2 <- Scenario_liste[, c("eir", "Access", "Dosage", "Indicator")]
  Scenario_liste_2 <- unique(Scenario_liste_2)

  # Calculate the number of mutation that have the same selection coefficient that emerge and establishment across each seed
  Scenario_liste_2$Number_mutation <- 0
  Scenario_liste_2$Number_mutation_2 <- 0
  Scenario_liste_2$esthablish <- 0

  for (this_s in Scenario_liste$Indicator) {
    Scenario_liste_2$Number_mutation[Scenario_liste_2$Indicator == this_s] <- sum(Scenario_liste$Number_mutation[Scenario_liste$Indicator == this_s], na.rm = T)
    Scenario_liste_2$Number_mutation_Celling[Scenario_liste_2$Indicator == this_s] <- sum(ceiling(Scenario_liste$Number_mutation[Scenario_liste$Indicator == this_s]), na.rm = T)
    Scenario_liste_2$esthablish[Scenario_liste_2$Indicator == this_s] <- sum(Scenario_liste$esthablisment[Scenario_liste$Indicator == this_s], na.rm = T)
  }

  # Estimate the probability of establishment
  Scenario_liste_2$Pe <- Scenario_liste_2$esthablish / Scenario_liste_2$Number_mutation
  Scenario_liste_2$Pe_celling  <- Scenario_liste_2$esthablish / Scenario_liste_2$Number_mutation_Celling
  
  # Save the files
  summary_path <- paste0(pm$pth$summarised, "Esthablishment.txt")
  write.table(Scenario_liste_2,
    file = summary_path, sep = "\t",
    col.names = TRUE, row.names = FALSE, quote = FALSE)
}

# ------------------------------------------------------------------------------------------------------------------
# Function to estimate number of mutation that emerge in the simulation by using the formula from Ian Hasting (2019).
# ------------------------------------------------------------------------------------------------------------------
Number_mutation <- function(Output_data_2) {

  # adjust the time step to year
  Output_data_2$Survey <- Output_data_2$Survey / Number_survey_years
  if (pm$opts$drug == "long") {
    Output_data_2 <- Output_data_2[Output_data_2$Survey >= 200, ]
    Output_data_2$Survey <- Output_data_2$Survey - 200
  } else {
    Output_data_2 <- Output_data_2[Output_data_2$Survey >= 250, ]
    Output_data_2$Survey <- Output_data_2$Survey - 250
  }

  # define the time at which we import mutation (mutation emerge)
  first_time <- Output_data_2$Survey[1]

  # control check to make sur all resistant genotype where eliminated before we start to import mutation
  if (Output_data_2$nInfectByGenotype_2[Output_data_2$Survey == first_time] >= 1) {
    Number_mut <- NA
  } else {

   # Estimate the number of time that a mutation emerge in the simulation based on importation rate and time of esthablisment
   Time_importation <- Output_data_2$Survey[Output_data_2$nInfectByGenotype_2 == 0]
   Time_esthablishment <- Time_importation[length(Time_importation)]
   Number_mut <- Scenario_liste$Importation[Scenario_liste$scenario_name == name_scenario] / 2 * 10 * Time_esthablishment
  }

  # return the number of mutation that emerged
  return(c(Number_mut))
}

# ----------------------------------------------------------------------------------------------------------------------------
# Function to estimate number of mutation that emerge in the simualtion by counting the number of time that a mutation emerge.
# ----------------------------------------------------------------------------------------------------------------------------
Number_mutation_2 <- function(Output_data_2) {

  # Adjust the timestep to year
  Output_data_2$Survey <- Output_data_2$Survey / Number_survey_years
  if (pm$opts$drug == "long") {
    Output_data_2 <- Output_data_2[Output_data_2$Survey >= 200, ]
    Output_data_2$Survey <- Output_data_2$Survey - 200
  } else {
    Output_data_2 <- Output_data_2[Output_data_2$Survey >= 250, ]
    Output_data_2$Survey <- Output_data_2$Survey - 250
  }

  # Define the time at which we start to import mutation
  first_time <- Output_data_2$Survey[1]

  # Control check to make sur all resistant genotype where eliminated before we start to import mutation
  if (Output_data_2$nInfectByGenotype_2[Output_data_2$Survey == first_time] >= 1) {
    Number_mut_2 <- NA
  } else {

    # Count the number of time that a mutation emerge in the simulation
    Data_number_resistant <- Output_data_2$nInfectByGenotype_2
    Number_mut_2 <- 1
    for (i in 2:length(Data_number_resistant)) {
      if (Data_number_resistant[i] == 0 & Data_number_resistant[i - 1] > 0) {
        Number_mut_2 <- Number_mut_2 + 1
      }
    }
  }

  return(c(Number_mut_2))
}

# ------------------------------------------------------------------
# Function to estimate the selection coefficient in perinal setting.
# ------------------------------------------------------------------
Esthablishment <- function(Output_data_2) {
  
  # Change the time in years
  Output_data_2$Survey <- Output_data_2$Survey / Number_survey_years
  
  # Delet data before we start to import mutation
  if (pm$opts$drug == "long") {
    Output_data_2 <- Output_data_2[Output_data_2$Survey >= 200, ]
    Output_data_2$Survey <- Output_data_2$Survey - 200
  } else {
    Output_data_2 <- Output_data_2[Output_data_2$Survey >= 250, ]
    Output_data_2$Survey <- Output_data_2$Survey - 250
  }

  # Define the tim zero at which we start to import mutation 
  first_time <- Output_data_2$Survey[1]
  
  # Control that there were no resistant genotype when we start to import mutation
  if (Output_data_2$nInfectByGenotype_2[Output_data_2$Survey == first_time] >= 1) {
    esthablish <- NA
  } else {
    
  # Estimate if one mutaiton establishe in this simulation
    if (Output_data_2$nInfectByGenotype_2[length(Output_data_2$nInfectByGenotype_2)] / Output_data_2$nInfect[length(Output_data_2$nInfectByGenotype_2)] >= 0.5) {
      esthablish <- 1
    } else {
      esthablish <- 0
    }
  }

  return(esthablish)
}

# ------------------------------------------------------------------
# Function to estimate the selection coefficient in perinal setting.
# ------------------------------------------------------------------
spread <- function(Output_data) {

  # Define the timestep
  time_step <- 30
  Number_survey_generation <- 60 / time_step # 1 generation equal 60 days
  Number_survey_years <- 365 / time_step

  # Define the initial time of the regression
  Time_regression <- Output_data$Survey[Output_data$nInfectByGenotype_2 / (Output_data$nInfectByGenotype_2 + Output_data$nInfectByGenotype) >= 0.5]
  Time_regression <- Time_regression[Time_regression >= 200 * Number_survey_years]
  time_start <- Time_regression[1]

  # Define the end time of the regression (2 year later)
  time_end <- time_start + 2 * Number_survey_years

  # Select the data within this boundary
  time_spread <- Output_data$Survey[Output_data$Survey > time_start & Output_data$Survey < time_end] / Number_survey_years
  Output_data$Inoculation_R <- Output_data$innoculationsPerAgeGroup_2 / (Output_data$innoculationsPerAgeGroup_2 + Output_data$innoculationsPerAgeGroup)
  Measurment_R <- Output_data$Inoculation_R[Output_data$Survey > time_start & Output_data$Survey < time_end]

  # Check that we have not to low frequency
  Measurment_R <- Measurment_R[Measurment_R <= 0.9 & Measurment_R >= 0.3]
  time_spread <- time_spread[1:length(Measurment_R)]

  # Estimate the selection coeffcient
  if (length(Measurment_R) >= 2 & is.na(sum(Measurment_R)) == F) {
    slope <- lm(log(Measurment_R / (1 - Measurment_R)) ~ time_spread)$coefficients[2]
    slope <- slope / 6
  } else {
    slope <- NA
  }

  return(slope)
}
