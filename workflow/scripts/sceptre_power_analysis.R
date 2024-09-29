# Script to create a sceptre object based off some simulated counts for downstream power simulations

### CREATE DEBUG FILES =======================================================

# Saving image for debugging
# save.image("RDA_objects/sceptre_power_analysis.rda")
# message("Saved Image")
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
  library(SingleCellExperiment)
  library(sceptre)
  source(file.path(snakemake@scriptdir, "R_functions/differential_expression_fun.R"))
  source(file.path(snakemake@scriptdir, "R_functions/power_simulations_fun.R"))
})

# Download Sceptre
library(devtools)
message("Installing Sceptre")
devtools::install_github("katsevich-lab/sceptre")
message("Sceptre Installation Complete")

# Load inputs
simulated_sceptre_object <- readRDS(snakemake@input$simulated_sceptre_object)
simulated_sce_disp <- readRDS(snakemake@input$simulated_sce_disp)
discovery_pairs_split <- read_tsv(snakemake@input$discovery_pairs_split, col_names = c("grna_group", "response_id"))
grna_target_data_frame <- read_tsv(snakemake@input$grna_target_data_frame)
response_matrix <- readRDS(snakemake@input$raw_counts)


# Load params
effect_size <- 1 - as.numeric(snakemake@params$effect_size)
reps <- snakemake@params$reps
guide_sd <- 0.13


### PRECOMPUTATIONS TO RUN SIMULATIONS ======================================

# Get all the perts to be run
perts <- unique(discovery_pairs_split$grna_group)

# Create an empty results data frame
discovery_results <- data.frame()


### RUN POWER SIMULATIONS ===================================================

for (pert in perts) {
 
  # Initialize the pert object with the given pert and sce
  pert_object <- pert_input(pert, sce = simulated_sce_disp, pert_level = "cre_perts")
  
  # Get all the guides that target the current `pert`
  pert_guides <- grna_target_data_frame %>%
    filter(grna_target == pert) %>%
    pull(grna_id)
  
  # Get perturbation status and gRNA perturbations for all cells
  pert_status <- colData(pert_object)$pert
  grna_perts <- assay(altExp(pert_object, "grna_perts"), "perts")
  # Convert to a sparse matrix, so the sampling function works in `create_guide_pert_status`
  grna_perts <- as(grna_perts, "CsparseMatrix")
  grna_pert_status <- create_guide_pert_status(pert_status, grna_perts = grna_perts, pert_guides = pert_guides)
  
  # Create effect size matrix (sampled from negative binomial distribution around effect_size or 1)
  effect_sizes <- structure(rep(effect_size, nrow(pert_object)), names = rownames(pert_object))
  
  # Loop through each rep
  for (rep in seq(reps)) {
    
    # Create and center effect size matrices
    es_mat <- create_effect_size_matrix(grna_pert_status, pert_guides = pert_guides,
                                        gene_effect_sizes = effect_sizes, guide_sd = guide_sd)
    es_mat <- center_effect_size_matrix(es_mat, pert_status = pert_status, gene_effect_sizes = effect_sizes)
    es_mat_use <- es_mat[, colnames(assay(pert_object, "counts"))]

    # Simulate Counts
    message("Simulating Counts")
    sim_counts <- sim_tapseq_sce(pert_object, effect_size_mat = es_mat_use)
    simulated_response_matrix <- as(assay(sim_counts, "counts"), "RsparseMatrix")


    # Save a temp sceptre object
    sceptre_object_use <- simulated_sceptre_object
    
    # Assign the simulated response matrix to the sceptre object
    sceptre_object_use@response_matrix <- list(simulated_response_matrix)
    # Subset and assign the discovery pairs relevant to the current perturbation to the sceptre object
    sceptre_object_use@discovery_pairs <- discovery_pairs_split[discovery_pairs_split$grna_group == pert,]
    
    # Run the discovery analysis
    message("Running discovery analysis")
    sceptre_object_use <- run_discovery_analysis(
      sceptre_object = sceptre_object_use,
      parallel = FALSE
    )

    # Get the discovery analysis results
    message("Returning discovery results")
    discovery_result <- get_result(
      sceptre_object = sceptre_object_use,
      analysis = "run_discovery_analysis"
    )

    # Add the number of perturbed cells and the rep to each pair
    n_pert_cells <- length(sceptre_object_use@grna_assignments$grna_group_idxs[[pert]])
    discovery_result$num_pert_cells <- n_pert_cells
    discovery_result$rep <- rep
    
    # Save the results
    discovery_results <- data.frame(rbind(discovery_results, discovery_result))

  }
}


### RUN POST COMPUTATIONS ===================================================

message("Processing output.")

# Construct a combined data frame
combined_data <- data.frame(
  response_id = rownames(simulated_sce_disp),
  disp_outlier_deseq2 = rowData(simulated_sce_disp)[, "disp_outlier_deseq2"],
  dispersion = rowData(simulated_sce_disp)[, "dispersion"],
  average_expression_all_cells = rowMeans(response_matrix),
  stringsAsFactors = FALSE
)

# Merge the combined data with discovery_results
discovery_results <- left_join(discovery_results, combined_data, by = "response_id")


### SAVE OUTPUT =============================================================

# save simulation output
message("Saving output to file.")
write_tsv(discovery_results, file = snakemake@output[[1]])


# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)