# Script to create a sceptre object based off some simulated counts for downstream power simulations

### CREATE DEBUG FILES =======================================================

# Saving image for debugging
save.image("RDA_objects/combine_sceptre_power_analysis.rda")
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Opening log file to collect all messages, warnings and errors
message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

### LOADING FILES ============================================================

# Load in the necessary packages
message("Loading packages")
suppressPackageStartupMessages({
  library(tidyverse)
})

### COMBINE SIMULATION OUTPUTS ==============================================

# Combine each file into one tsv and save it
message("Combining power analysis outputs...")

# Initialize an empty list to store data frames
data_list <- list()

# Loop over input files
message("Gathering inputs")
for (file_path in snakemake@input) {
  # Read the current file
  current_data <- read_tsv(file_path, col_types = cols())
  
  # Append the data frame to the list
  data_list[[length(data_list) + 1]] <- current_data
}

# Combine all data frames into one
message("Combining all data frames into one")
combined_data <- do.call(rbind, data_list)

### SAVE COMBINE OUTPUTS ===================================================

# Write the combined data frame to the output file
write_tsv(combined_data, snakemake@output[[1]])
message("Power analysis outputs combined successfully.")

# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)
