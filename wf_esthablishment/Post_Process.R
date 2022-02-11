######################################################################################
#     Post process of the data                                                       #
#                                                                                    #
# Task :Reorganize the output of OpenMalaria in a table                              #
#       that display all the measurement at each time step                           #
#                                                                                    #
# Input: Outpout file form Openmalaria                                               #
#                                                                                    #
# Output: Table of each measurement values at each time step                         #
#                                                                                    #
# author: thiery.masserey@swisstph.ch                                                #
######################################################################################

Postprocess <- function(pm) {
  # Message in the console
  message("  - Post Processing ")
  
  # Define list of output files
  Output <- list.files(pm$pth$sim_out) # Make a list of all the output file
  
  # Initiate progress bar
  pb = txtProgressBar(
    min = 0,
    max = length(Output),
    initial = 0,
    width = 100,
    style = 3)
  
  # Start the loop that will do the post process for each output file
  for (i in 1:length(Output)) {
    
    # define file name, file path, and download it i<-1
    Output_data_name <- Output[i]  # select the file name
    Output_data_file_path <- file.path(pm$pth$sim_out, Output_data_name) # define the path to the file
    Output_data <- read.table(Output_data_file_path) # open the file
    
    # Extract the survey Time in a new dataframe (Output_data_2)
    a <- as.data.frame(table(Output_data$V1))
    Survey <- as.numeric(a$Var1)
    Output_data_2 <- as.data.frame(Survey)
    
    # Extract the different Measure and age group and define the meaning of the measurement
    Indicator <- Output_data$V3[Output_data$V1 == 1] # indicators
    Age_groupe <- Output_data$V2[Output_data$V1 == 1] # age group
    # Indicator meaning if look at ACT
    if (pm$opts$drug == "ACT") {
      Indicators_meaning <-
        c("nHost",            #                 0:  Total number of humans
          "nInfect",          #                 1:  The number of human hosts with an infection
          "nExpectd",          #                2:  The expected number of infected hosts
          "nPatent",          #                 3:  The number of human hosts whose total (blood-stage) parasite density is above the detection threshold
          "nTransmit",          #               7:  Infectiousness of human population to mosquitoes: sum(p(transmit_i)) across humans i, weighted by availability to mosquitoes.
          "nTreatments1",          #            11: Number of treatments (1st line)
          "nUncomp",          #                 14: Number of episodes uncomplicated1
          "innoculationsPerAgeGroup",         # 30: The total number of inoculations per genotype (genotype 1)
          "innoculationsPerAgeGroup_2",       # 30: The total number of inoculations per genotype (genotype 2)
          "Vector_Nv",          #               32: Host seeking mosquito population size at this time step, species 1
          "Vector_Nv_2",          #             32: Host seeking mosquito population size at this time step, species 2
          "inputEIR",          #                35: Simulated EIR
          "simulatedEIR",          #            36: Input EIR
          "nMDAs",          #                   52: Number of drug doses given via mass deployment (MDA or screen&treat)
          "nMassScreenings",          #         55: Report number of screenings used by MDA/MSAT
          "nTreatDiagnostics",          #       64: Number of diagnostic tests performed (if in the health system description, useDiagnosticUC="true").
          "nInfectByGenotype",          #       69: The number of human hosts with an infection (patent or not), for each genotype (genotype 1)
          "nInfectByGenotype_2",          #     69: The number of human hosts with an infection (patent or not), for each genotype (genotype 2)
          "nPatentByGenotype",          #       70: The number of human hosts whose total (blood-stage) parasite density, for each genotype, is above the detection threshold (genotype 1)
          "nPatentByGenotype_2",          #     70: The number of human hosts whose total (blood-stage) parasite density, for each genotype, is above the detection threshold (genotype 2)
          "nHostDrugConcNonZero",          #    72: For each drug type in the pharmacology section of the XML, report the number of humans with non-zero concentration
          "nHostDrugConcNonZero_2",
          "nHostDrugConcNonZero_3",
          "nHostDrugConcNonZero_4",
          "nHostDrugConcNonZero_5",
          "nHostDrugConcNonZero_6",
          "sumLogDrugConcNonZero",          #   73: For each drug type in the pharmacology section of the XML, report the sum of the natural logarithm of the drug concentration in hosts with non-zero concentration.
          "sumLogDrugConcNonZero_2",
          "sumLogDrugConcNonZero_3",
          "sumLogDrugConcNonZero_4",
          "sumLogDrugConcNonZero_5",
          "sumLogDrugConcNonZero_6")
      
    } else{
      # Indicator meaning if look at short or long acting drug
        Indicators_meaning <-
        c("nHost",            #                 0:  Total number of humans
          "nInfect",          #                 1:  The number of human hosts with an infection
          "nExpectd",          #                2:  The expected number of infected hosts
          "nPatent",          #                 3:  The number of human hosts whose total (blood-stage) parasite density is above the detection threshold
          "nTransmit",          #               7:  Infectiousness of human population to mosquitoes: sum(p(transmit_i)) across humans i, weighted by availability to mosquitoes.
          "nTreatments1",          #            11: Number of treatments (1st line)
          "nUncomp",          #                 14: Number of episodes uncomplicated1
          "innoculationsPerAgeGroup",         # 30: The total number of inoculations per genotype (genotype 1)
          "innoculationsPerAgeGroup_2",       # 30: The total number of inoculations per genotype (genotype 2)
          "Vector_Nv",          #               32: Host seeking mosquito population size at this time step, species 1
          "Vector_Nv_2",          #             32: Host seeking mosquito population size at this time step, species 2
          "inputEIR",          #                35: Simulated EIR
          "simulatedEIR",          #            36: Input EIR
          "nMDAs",          #                   52: Number of drug doses given via mass deployment (MDA or screen&treat)
          "nMassScreenings",          #         55: Report number of screenings used by MDA/MSAT
          "nTreatDiagnostics",          #       64: Number of diagnostic tests performed (if in the health system description, useDiagnosticUC="true").
          "nInfectByGenotype",          #       69: The number of human hosts with an infection (patent or not), for each genotype (genotype 1)
          "nInfectByGenotype_2",          #     69: The number of human hosts with an infection (patent or not), for each genotype (genotype 2)
          "nPatentByGenotype",          #       70: The number of human hosts whose total (blood-stage) parasite density, for each genotype, is above the detection threshold (genotype 1)
          "nPatentByGenotype_2",          #     70: The number of human hosts whose total (blood-stage) parasite density, for each genotype, is above the detection threshold (genotype 2)
          "nHostDrugConcNonZero",          #    72: For each drug type in the pharmacology section of the XML, report the number of humans with non-zero concentration
          "nHostDrugConcNonZero_2",
          "nHostDrugConcNonZero_3",
          "nHostDrugConcNonZero_4",
          "nHostDrugConcNonZero_5",
          "sumLogDrugConcNonZero",          #   73: For each drug type in the pharmacology section of the XML, report the sum of the natural logarithm of the drug concentration in hosts with non-zero concentration.
          "sumLogDrugConcNonZero_2",
          "sumLogDrugConcNonZero_3",
          "sumLogDrugConcNonZero_4",
          "sumLogDrugConcNonZero_5")
    }
    
    # Extract the information for each measurement, and save it in the new dataframe (Output_data_2)
    for (j in 1:(length(Indicator))) {
      Output_data_2[j + 1] <- Output_data$V4[Output_data$V3 == Indicator[j] & Output_data$V2 == Age_groupe[j]] # put in the j +1 colone (due that we have one colone of the time survey), the value of the jth indicator
      colnames(Output_data_2)[j + 1] <- Indicators_meaning[j] # give the name of the indicators
    }
    
    # Save the new dataset in the post process folder
    Postprocess_data_name <- paste0('PostProcess_', Output_data_name)
    Postprocess_data_path <- file.path(pm$pth$processed, Postprocess_data_name)
    write.table(Output_data_2, file = Postprocess_data_path, sep = ";")
    data_processed_5day <- read.table(file = Postprocess_data_path, sep = ";")
    
    # Update progress bar
    setTxtProgressBar(pb, i)
    
  } # close the loop
  
  close(pb) # close the progress bar
  
} # close the function
