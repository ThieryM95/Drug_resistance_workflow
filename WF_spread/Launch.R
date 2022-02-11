#######################################################################################
# Workflow to estimate the impact of parameter on the selection coefficient           #
# via global sensitivity                                                              #
#                                                                                     #
# Input: Need a basefile, defined the parameter of interest and their range           #
# Aim: Assess the impact of each parameter on the rate of spread of resistant         #
#      genotype                                                                       #
# Output: Sobol indices of each factor on the rate of spread, and direction           #
#         of the impact of each factor on the rate of spread                          #
#                                                                                     #
# Run on OpenMalaria V40.1                                                            #
# author: thiery.masserey@swisstph.ch adapted from andrewjames.shattock@swisstph.ch   #
#######################################################################################

# Clear global environment
rm(list = ls())

# Set the working directory
setwd("/scicore/home/penny/masthi00/WF_spread")

# download functions in the environment
source("myRfunctions.R")
source("directories.R")
source("openmalaria_setup_seed.R")
source("openmalaria_run.R")
source("Post_Process.R")
source("summary_results.R")
library("tgp")
library("styler")
source("Option.R")
source("Generate_parameter_table.R")
source("GP.R")
source("adaptative_sample.R")
source("sensitivity_analysis.R")
source("Post_process_GP.R")
source("Post_process_sensitivity_analysis.R")

# Tidy up
# Clear console
if (interactive()) clc() # see myRfunctions.R

## Close figures
if (interactive()) clf() # see myRfunctions.R

#########################
# Start of the workflow #
#########################

# 1) Creat directory, and define options (define which drug archetype: A, B or A+B):
pm <- set_options() # see Option.R

# 2) Generate latin hypercube sample and list of scenario -
generate_param_table(pm, param_table = as_param_table) # see Generate_parameter_table.R

# 3) Generate seed paterns
n_jobs_unique <- simulation_sed_patterns(pm) # see openmalaria_setup_seed.R
message("   ~ Total number of scenarios: ", thou_sep(n_jobs_unique), " (with seed)")

# 4) Generate scenario xml files from the base xml
generate_xml_files(
  pm,
  n_jobs = n_jobs_unique,
  file_path = pm$pth$xml_base,
  sed_path  = pm$pth$sim_sed,
  xml_path  = pm$pth$sim_xml) # see openmalaria_setup_seed.R

# 5) Run simulation
run_model(
  pm,
  n_jobs_unique,
  xml_path = pm$pth$sim_xml,
  sim_path = pm$pth$sim_out) # see openmalaria_run.R

# 6) Post-processing
Postprocess(pm) # see Post_Process.R

# 7) Summarized the results
SummaryResults(pm) # see summary_results.R

# 8) fit the Gaussian process
Results_gp_1 <- run_gp(pm) # see  GP.R

# 9) Start adaptive sampling
Results_gp <- run_adaptive_sampling(pm) # see adaptative_sample.R

# 10) Post process to visualize the results of the GP in a table format
Precision_final <- Post_process_GP(Results_gp) # see Post_process_GP.R

# 11)  Sensitivity analysis
Results_SA <- sensitivity_analysis(pm) # see sensitivity_analysis.R

# 12) Post process the result of the sensitivity analysis in a table format

# Function to transform the sobol indies into a table format
data <- Post_process_sensitivity(Results_SA) # see Post_process_sensitivity_analysis.R

# Function to transform the direction of the effect of parameter on the selection coeffcient in a table format
Quantil_final_final <- Post_process_sensitivity_2(Results_SA) # see Post_process_sensitivity_analysis.R

# See the Visualize_results folder to Visualize the results