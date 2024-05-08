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


### SAVE OUTPUT ============================================================

# save output sce to file
write_tsv(power, file = snakemake@output$power_analysis_results)

# close log file connection
sink()
sink(type = "message")
close(log)



