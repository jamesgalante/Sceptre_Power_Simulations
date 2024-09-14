## Calculate the power for each pair from the power simulations output

# Saving image for debugging
save.image("RDA_objects/compute_power_from_simulations.rda")
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# opening log file to collect all messages, warnings and errors
message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


### LOADING INPUTS ==========================================================

# required packages and functions
suppressPackageStartupMessages({
  library(tidyverse)
})

# Load input file
power_sim <- read_tsv(snakemake@input$power_analysis_output)


### CALCULATE POWER WITH ONE OF TWO METHODS
if (is.null(snakemake@params$positive_proportion)) {
### CALCULATING POWER ======================================================

# Calculate power
power <- power_sim %>%
  group_by(rep) %>%
  mutate(pval_adj = p.adjust(p_value, method = "BH")) %>%
  mutate(pval_adj = replace_na(pval_adj, 1)) %>%
  group_by(grna_target, response_id, disp_outlier_deseq2, average_expression_all_cells) %>%
  summarize(
    mean_log2_fold_change = mean(log_2_fold_change, na.rm = TRUE),  # Calculate mean ignoring NAs
    mean_pert_cells = mean(num_pert_cells, na.rm = TRUE),  # Calculate mean perturbed cells, ignoring NAs
    power = mean(pval_adj < snakemake@params$pval_adj_thresh & (replace_na(log_2_fold_change, 1) < 0), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(power))

# Add the effect size (NOTE: if adding ability to have multiple effect sizes in one run, must change this to a wildcard)
power$effect_size <- as.numeric(snakemake@params$effect_size)

} else {
### CALCULATING POWER WITH CORRECT POSITIVES PROPORTION ====================

# Initialize an empty list to store results for each rep
power_results <- list()

# Calculate the number of null p-values to sample (constant for all reps)
num_rows_first_rep <- nrow(power_sim[power_sim$rep == 1, ])
positive_proportion <- as.numeric(snakemake@params$positive_proportion)
n_null <- round(num_rows_first_rep * ((1 - positive_proportion) / positive_proportion))

# Loop over each rep, calculating adjusted p values with positive proportion in mind
for (i in unique(power_sim$rep)) {
  # Grab data for the current rep
  rep_data <- power_sim %>% filter(rep == i)
  
  # Sample null p-values from uniform(0, 1)
  null_pvals <- runif(n_null, min = 0, max = 1)
  
  # Combine real and null p-values
  combined_pvals <- c(rep_data$p_value, null_pvals)
  
  # Perform Benjamini-Hochberg correction on the combined p-values
  pval_adj_combined <- p.adjust(combined_pvals, method = "BH")
  
  # Keep only the first `num_rows_first_rep` adjusted p-values (real ones)
  rep_data$pval_adj <- pval_adj_combined[1:num_rows_first_rep]
  
  # Replace any NA adjusted p-values with 1
  rep_data$pval_adj <- replace_na(rep_data$pval_adj, 1)
  
  # Calculate mean log2 fold change and other metrics
  rep_summary <- rep_data %>%
    group_by(grna_target, response_id, disp_outlier_deseq2, average_expression_all_cells) %>%
    summarize(
      mean_log2_fold_change = mean(log_2_fold_change, na.rm = TRUE),  # Calculate mean ignoring NAs
      mean_pert_cells = mean(num_pert_cells, na.rm = TRUE),  # Calculate mean perturbed cells, ignoring NAs
      power = mean(pval_adj < snakemake@params$pval_adj_thresh & (replace_na(log_2_fold_change, 1) < 0), na.rm = TRUE),
      .groups = "drop"
    )
  
  # Append the result to the list
  power_results[[i]] <- rep_summary
}

# Combine the results into one data frame and add effect size
power <- bind_rows(power_results)
power$effect_size <- as.numeric(snakemake@params$effect_size)

# Sort the results by power
power <- power %>%
  arrange(desc(power))

}
### SAVE OUTPUT ============================================================

# save output sce to file
write_tsv(power, file = snakemake@output$power_analysis_results)

# close log file connection
sink()
sink(type = "message")
close(log)



