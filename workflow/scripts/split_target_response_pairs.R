
# Saving image for debugging
save.image("RDA_objects/split_target_response_pairs.rda")
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")


# opening log file to collect all messages, warnings and errors
message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

### LOAD DATA ===============================================================

# Load packages
message("Loading packages")
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tidyverse)
})

# Retrieve parameters from Snakemake
message("Loading inputs")
simulated_sce_disp <- readRDS(snakemake@input$simulated_sce_disp)
batches <- as.numeric(snakemake@params$batches)
output_files <- snakemake@output$splits


### CREATE DISCOVERY PAIRS FILE =============================================

message("Creating Discovery Pairs file")
grna_targets <- rownames(altExp(simulated_sce_disp, "cre_perts"))
response_ids <- rownames(simulated_sce_disp)

# Remove NA dispersion responses (resulting from rowSum(gene) < 2) as these counts won't be simulated
response_ids <- response_ids[!is.na(rowData(simulated_sce_disp)[, "dispersion"])]

# Create and save discovery_pairs
discovery_pairs <- expand.grid(grna_target = grna_targets, response_id = response_ids)
write_tsv(discovery_pairs, snakemake@output$discovery_pairs)


### CREATE GRNA_TARGET_DATA_FRAME ===========================================

message("Creating grna_target_data_frame")
grna_target_data_frame <- data.frame(
  grna_id = rowData(altExp(simulated_sce_disp, "grna_perts"))$name,
  grna_target = rowData(altExp(simulated_sce_disp, "grna_perts"))$target_name
)
write_tsv(grna_target_data_frame, snakemake@output$grna_target_data_frame)


### SPLIT DISCOVERY_PAIRS FILE ==============================================

# Calculate the target number of unique grna_group for each split
message("Precomputations for splitting discovery_pairs")
total_unique_groups <- n_distinct(discovery_pairs$grna_target)
target_per_split <- ceiling(total_unique_groups / batches)

# Make sure that splitting the unique groups will work given the number of batches requested
if (total_unique_groups < batches) {
  stop("The number of batches for parallelization exceeds the number of unique grna_targets.")
}

# Initialize splits with empty data frames
splits <- vector("list", batches)
names(splits) <- paste0("split", seq_len(batches))
for (i in seq_len(batches)) {
  splits[[i]] <- data.frame(grna_target = character(0), response_id = character(0))
}

# Distribute grna_target to splits trying to even out the number of unique values
message("Distributing discovery_pairs")
for (i in seq_len(total_unique_groups)) {
  split_counts <- sapply(splits, function(x) n_distinct(x$grna_target))
  split_with_least <- which.min(split_counts)
  
  # Get the rows for the current grna_target
  current_rows <- discovery_pairs %>% 
    filter(grna_target == unique(discovery_pairs$grna_target)[i])
  
  # Add the current rows to the appropriate split
  splits[[split_with_least]] <- bind_rows(splits[[split_with_least]], current_rows)
}

# Write each split to the corresponding output file
message("Saving splits to files")
for (i in seq_along(splits)) {
  write_tsv(splits[[i]], file=output_files[[i]], col_names = FALSE)
}


# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)
