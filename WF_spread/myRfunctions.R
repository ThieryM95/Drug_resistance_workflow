####################################################
# A series of helpful R functions.                ##
#                                                 ##
# author: andrewjames.shattock@swisstph.ch        ##
# adapted: thiery.masserey@swisstph.ch            ##
####################################################

# load package
require(rlist)
require(stringr)
require(fastmatch)
require(zoo)

# ---------------------------------------------------------
# Assert that a condition is true. Throw error if not.
# ---------------------------------------------------------
assert <- function(condition, errMsg) {
  if (!condition) {
    stop(errMsg)
  }
}

# ---------------------------------------------------------
# Clear the console.
# ---------------------------------------------------------
clc <- function() {
  cat("\014")
}

# ---------------------------------------------------------
# Clear all figures.
# ---------------------------------------------------------
clf <- function() {
  graphics.off()
}

# ---------------------------------------------------------
# Clear environment (and also console and figures by default).
# ---------------------------------------------------------
clear <- function(clearConsole = TRUE, clearFigs = TRUE) {
  if (clearConsole == TRUE) {
    clc()
  }
  if (clearFigs == TRUE) {
    clf()
  }
  rm(list = ls(envir = globalenv()), envir = globalenv())
}

# ---------------------------------------------------------
# Is an element conatined within a vector? Returns boolean.
# ---------------------------------------------------------
contains <- function(vec, el) {
  bool <- is.element(el, vec) # Simply for readability
  return(bool)
}

# ---------------------------------------------------------
# Create a log file (see function wait_for_jobs)
# ---------------------------------------------------------
create_bash_log <- function(log_file) {
  if (file.exists(log_file)) file.remove(log_file)
  Sys.sleep(0.1)
  file.create(log_file)
  Sys.sleep(0.1)
}

# ---------------------------------------------------------
# Create dataframe with custom column headings in one line.
# ---------------------------------------------------------
dataFrame <- function(x, colNames = NULL, rowNames = NULL, ...) {
  df <- data.frame(x, row.names = rowNames, ...) # Create dataframe
  if (!is.null(colNames)) {
    colnames(df) <- colNames
  } # Rename column headers if need be
  return(df)
}

# ---------------------------------------------------------
# Drop one or more columns from a dataframe - improved readability.
# ---------------------------------------------------------
df_drop <- function(df, drop_cols) {
  drop_idx <- (names(df) %in% drop_cols)
  df <- df[, !drop_idx]
  return(df)
}

# ---------------------------------------------------------
# Convert dataframe to a list using a column as list fields.
# ---------------------------------------------------------
df2list <- function(df, col_idx, append_list = NULL) {
  if (is.null(append_list)) append_list <- list() # Initiate list if needed
  for (i in 1:nrow(df)) {
    vec <- unlist(unname(df[i, -col_idx])) # Convert to vector
    vec <- na.trim(vec, sides = "right") # Remove trailing NAs
    append_list[[df[i, col_idx]]] <- vec # Append vector to list
  }
  return(append_list)
}

# ---------------------------------------------------------
# Create figure window and display plot.
# ---------------------------------------------------------
figure <- function(figHandle) {
  operatingSystem <- Sys.info()["sysname"]
  if (operatingSystem == "Linux") plot(figHandle)
  if (operatingSystem == "Windows") {
    windows() # Other techniques available
    print(figHandle)
  }
}

# ---------------------------------------------------------
# Platform specific file separator - for readability.
# ---------------------------------------------------------
file_sep <- function() {
  platform_file_sep <- .Platform$file.sep
  return(platform_file_sep)
}

# ---------------------------------------------------------
# Display blank line in console.
# ---------------------------------------------------------
leaveLine <- function(nLines = 1) {
  for (i in 1:nLines) {
    message("")
  }
}

# ---------------------------------------------------------
# Convert a field of a list into a vector.
# ---------------------------------------------------------
list2vec <- function(l, idx = NA) {
  vals <- unlist(lapply(l, "[[", idx))
  return(vals)
}

# ---------------------------------------------------------
# Number of running and pending jobs on the cluster
# ---------------------------------------------------------
n_slurm_jobs <- function(user = "shatto0000") {

  # Base sq command for user
  sq <- paste("squeue -u", user)

  # Concatenate full commands
  slurm_running <- paste(sq, "-t running | wc -l")
  slurm_pending <- paste(sq, "-t pending | wc -l")

  # System call to determine number of slurm processes
  n_running <- system(slurm_running, intern = TRUE)
  n_pending <- system(slurm_pending, intern = TRUE)

  # Convert to numeric and minus 1 to get number of jobs
  n_running <- as.numeric(n_running) - 1
  n_pending <- as.numeric(n_pending) - 1

  # Compile into list
  n_jobs <- list(running = n_running, pending = n_pending)

  return(n_jobs)
}

# ---------------------------------------------------------
# Sum whilst ignoring NAs - better readability.
# ---------------------------------------------------------
nasum <- function(x) {
  y <- sum(x, na.rm = TRUE)
  return(y)
}

# ---------------------------------------------------------
# Suppress output from a function call.
# ---------------------------------------------------------
quiet <- function(x) {
  sinkCon <- file("sink.txt")
  sink(sinkCon, type = "output")
  sink(sinkCon, type = "message")
  on.exit(sink(type = "output"))
  on.exit(sink(type = "message"), add = TRUE)
  on.exit(file.remove("sink.txt"), add = TRUE)
  invisible(force(x))
}

# ---------------------------------------------------------
# Concatenate strings with no seperator.
# ---------------------------------------------------------
strcat <- function(...) {
  str <- paste(..., sep = "")
  return(str)
}

# ---------------------------------------------------------
# Sumbit an array job to the SciCORE cluster.
# ---------------------------------------------------------
submit_cluster <- function(bash_file, n_jobs, ...) {

  # Create a new log file for the cluster jobs
  log_file <- file.path("scicore_log.txt")
  create_bash_log(log_file)

  # Construct slurm array command for running in parallel
  slurm_array <- paste0("--array=1-", n_jobs)

  # Concatenate system command
  sys_command <- paste(
    "sbatch", slurm_array,
    bash_file, log_file, ...
  )

  # Invoke this command
  system(sys_command)

  browser() # Need to rejig wait_for_jobs a little...

  # Wait for all cluster jobs to complete
  wait_for_jobs(pm, log_file, n_jobs)
}

# ---------------------------------------------------------
# Summing function for higher dimensional arrays.
# ---------------------------------------------------------
sumMat <- function(x, sumDim) {
  dims <- dim(x)
  nDims <- length(dims)
  if (nDims == 2) { # 2D matrices included for readability
    if (sumDim == 1) {
      y <- rowSums(x) # Simply apply rowSums
    } else if (sumDim == 2) {
      y <- colSums(x) # Simply apply colSums
    }
  }
  if (nDims > 2) { # Higher dimensional arrays are summed and squeezed
    applyDim <- setdiff(1:nDims, sumDim)
    y <- apply(x, applyDim, sum) # Use apply method
  }
  return(y)
}

# ---------------------------------------------------------
# Format a number with thousand mark separators.
# ---------------------------------------------------------
thou_sep <- function(val) {
  format_val <- format(val, big.mark = ",") # Automatically rounds
  return(format_val)
}

# ---------------------------------------------------------
# Repeat vector into a 2D matrix.
# ---------------------------------------------------------
vec2mat <- function(v, nRow = 1, nCol = 1) {
  m <- matrix(v, nrow = nRow, ncol = nCol * length(v), byrow = TRUE)
  return(m)
}

# ---------------------------------------------------------
# Waits until all jobs have finished (marked in a log file).
#
# Adapted from code written by Monica Golumbeanu.
# ---------------------------------------------------------
wait_for_jobs <- function(pm, log_file, n_lines, wait_time = 1) {

  # message(" - Waiting for batch job to complete...")

  # Assume we won't have any issues
  success_flag <- TRUE

  # Wait for log file to be created
  while (!file.exists(log_file)) {
    Sys.sleep(wait_time)
  }

  # Initiate a progress bar
  pb <- txtProgressBar(
    min = 0,
    max = n_lines,
    initial = 0,
    width = 100,
    style = 3
  )

  # Wait for at least one line to be written
  while (file.info(log_file)$size == 0) {
    Sys.sleep(wait_time)
  }

  # TODO: Best approach is to wait until n_jobs == 0
  #       Then check if all jobs have been successfully created

  # Wait for scenarios to be created
  while (nrow(read.table(log_file)) < n_lines) {

    # Update progress bar
    setTxtProgressBar(pb, nrow(read.table(log_file)))

    # Check number of running and pending jobs
    n_jobs <- n_slurm_jobs(user = pm$user)

    # If none are left, then we've had some failures
    if (sum(unlist(n_jobs)) == 0) {
      success_flag <- FALSE

      break # Break out of while loop
    }

    # Wait before testing again
    Sys.sleep(wait_time)
  }

  # Finalise progress bar
  setTxtProgressBar(pb, n_lines)
  close(pb)

  if (success_flag == TRUE) {
    # message(" - Batch job complete!")

    return(TRUE)
  } else {
    message("!! Batch job finished with errors !!")

    # Check which jobs have failed
    # failed_jobs = check_jobs(pm)

    # return(failed_jobs)
  }
}

# ---------------------------------------------------------
# Writes a line of text to file (handling file connection).
# ---------------------------------------------------------
write_lines <- function(txt, file_name) {
  file_conn <- file(file_name, blocking = TRUE)
  writeLines(txt, file_conn)
  close(file_conn)
}

#----------------------------------------------------
# Adjuste the EIR value (needed in version 40.1)
#----------------------------------------------------
# In the perineal setting
adjust_EIR <- function(eir, access) {

  # load the emulator
  trained_model <- readRDS(paste0(file.path("/scicore/home/penny/masthi00", "WF_spread", "EIR_gp_no_seasonality.RData")))

  # predict the input EIR to simulat the wanted EIR in a setting that has a specific treatment access
  predict_eir <- predict(x = as.matrix(data.frame(eir, access)), object = trained_model)
  predict_eir <- predict_eir$mean

  # return the adjusted eir
  return(predict_eir)
}

# In the seasonal setting
adjust_EIR_seasonality <- function(eir, access) {

  # load the emulator
  trained_model <- readRDS(paste0(file.path("/scicore/home/penny/masthi00", "WF_spread", "EIR_gp_seasonality.RData")))

  # Predict the input EIR to simulat the wanted EIR in a setting that has a specific treatment access
  predict_eir <- predict(x = as.matrix(data.frame(eir, access)), object = trained_model)
  predict_eir <- predict_eir$mean

  # return the adjusted eir
  return(predict_eir)
}

#---------------------------------------------------------------------------
# Estimate the maximum killing rate needed to model a specific fitness cost
#---------------------------------------------------------------------------
adjust_Fitness <- function(fitness) {
  
  # load the emulator
  trained_model <- readRDS(paste0(file.path("/scicore/home/penny/masthi00", "WF_spread", "Fitness_gp.RData")))
  
  # predict the maximum killing rate needed to model a specific fitness cost
  NEW_DATA <- NULL
  NEW_DATA$Scenario_liste.Indicator <- fitness
  predict_eir <- predict(trained_model, newdata = NEW_DATA)
  
  # return the maximum killing rate
  return(predict_eir)
}

#---------------------------------------------------
# Define parameter table name based on the iteration
#--------------------------------------------------- 
param_set_name <- function(sample_num = NA, all_samples = FALSE) {

  # Default name for all samples
  if (all_samples == TRUE) {
    name <- "parameter_table_all_samples"
  } else {

    # Set in one place so can be more easily changed
    name <- paste0("parameter_table_sample_", sample_num)
  }

  # return the parameter name
  return(name)
}
