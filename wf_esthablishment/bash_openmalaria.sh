#!/bin/bash
#SBATCH --job-name=run_OpenMalaria
#SBATCH --account=penny
#SBATCH --time=05:59:00  # 03:00:00
#SBATCH --qos=6hours     #30min 6hours
#SBATCH --mem=10G         # 4G 10G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scicore/home/penny/masthi00/wf_esthablishment/JOB_OUT/scicore_output%A%a.txt
#############################################################################################
# Script for submitting OpenMalaria simulations as jobs to the cluster.                     #
#                                                                                           #
#                                                                                           #
# Input : Xml file for each scenario                                                        #
#                                                                                           #
# Task: Run each scenario with OpenMalariaon on the cluster                                 #
#                                                                                           #
# Output: One scenario output file for each simulation                                      #                                                    #
#                                                                                           #
# Arguments:                                                                                #
#		 INPUT_DIR: Directory containing the scenario xml files                             #
#		  DEST_DIR: Directory where two output files will be created per OpenMalaria        #
#                   simulation                                                              #
#	  RESOURCE_DIR: Directory of OpenMalaria resources                                      #
#                                                                                           #
# Created by monica.golumbeanu@swisstph.ch  and andrewjames.shattock@swisstph.ch            #                                                      #
# Adapted by thiery.masserey@swisstph.ch                                                    #
#############################################################################################

# Load OpenMalaria module
module purge
ml  OpenMalaria/40.1-iomkl-2019.01

# Define input variables
INPUT_DIR=$1
DEST_DIR=$2
RESOURCE_DIR=$3
LOG_FILE=$4

# Change to folder with resource files
cd $RESOURCE_DIR

# Extract scenario xml files
SCENARIO_FILES=(${INPUT_DIR}*.xml)

# Select single scenario file based on task ID
ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)
SCENARIO_FILE=${SCENARIO_FILES[$ID]}

# Construct paths to xml file and output dirs
SCENARIO_NAME=$(basename "$SCENARIO_FILE" ".xml")
OUTPUT1=$DEST_DIR$SCENARIO_NAME"_out.txt"
OUTPUT2=$DEST_DIR$SCENARIO_NAME"_cts.txt"

# Check if the output file already exists
if [ -f "$OUTPUT1" ]; then
    echo "Scenario $SCENARIO_NAME has already been run - not re-running"
    
else 
    echo "Running simulation for $SCENARIO_NAME"

    # Run OpenMalaria for this scenario file
    openMalaria --scenario $SCENARIO_FILE --output $OUTPUT1 --ctsout $OUTPUT2
fi

# Write scenario name to log file
echo "OpenMalaria simulation $SCENARIO_NAME finished" >> $LOG_FILE

