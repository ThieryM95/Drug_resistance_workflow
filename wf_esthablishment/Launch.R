################################################################################
# Workflow to estimate the probability of establishment of mutation that       #
# have a specific selection coefficient                                        #                                                     #                        ##
#                                                                              #
# Task: estimate the relationship between selection coefficient and probability#
#       of establishment                                                       #
#                                                                              #
# Input: Need a base file, define the settings, and table of parameter value   #
#        and selection coefficient                                             #
#                                                                              #
# Output: Relation between probability of establishment and selection          #
#         coefficient                                                          #
#                                                                              #
# Run on OpenMalaria V40.1                                                     #
#                                                                              #
# author: Thiery.masserey@swisstph.ch                                          #
################################################################################

# Clear global environment
rm(list = ls())

# Set working directory to sourced file
setwd("/scicore/home/penny/masthi00/wf_esthablishment")

source("myRfunctions.R")
source("directories.R")
source("Option.R")
source("Generate_paramt_table.R")
source("openmalaria_setup_seed.R")      
source("openmalaria_run.R")
source("Post_Process.R")
source("summary_results.R")
library("styler")
library("tgp")

# Tidy up
# Clear console
if (interactive()) clc()

# Close figures
if (interactive()) clf()

#########################
# Start of the workflow #
#########################

# 1) Creat directory, and define options (define which drug archetype: A, B or A+B):
pm <- set_options(do_step = 4) # see Option.R

# 2) Generate parameter table by selecting combination of parameter that have a know selection coefficient
generate_param_table(pm, param_table = as_param_table) # see generate_parameter_table.R

# 3) Generate seed patterns
n_jobs_unique <- simulation_sed_patterns(pm) # see open_malaria_setup_seed.R
message("   ~ Total number of scenarios: ", thou_sep(n_jobs_unique), " (with seed)")

# 4) Generate scenario xml files from the base xml
generate_xml_files(pm,
  n_jobs = n_jobs_unique, 
  file_path = pm$pth$xml_base,
  sed_path = pm$pth$sim_sed,
  xml_path = pm$pth$sim_xml) # see open_malaria_setup_seed.R

# 5) Run simulation
run_model(pm, n_jobs_unique,
  xml_path = pm$pth$sim_xml,
  sim_path = pm$pth$sim_out) # see openmalaria_run.R

# 6) Post processing
Postprocess(pm) # see Post_process.R

# 7) Summarized the results
SummaryResults(pm) # see Summary_results.R

# See visualise_results folder to visualize the results