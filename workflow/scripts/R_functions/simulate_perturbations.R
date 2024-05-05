
## Functions to simulate enhancer perturbations by randomly picking cells perturbed for a specified
## number of enhancers

library(SingleCellExperiment)

# function to simulate enhancer perturbations for randomly picked cells and add as perturbation
# status matrix to the SCE object
simulate_perturbations <- function(sce, cells_per_pert, guides_per_pert) {
  
  # create randomly selected perturbations
  cells <- colnames(sce)
  pert_status <- lapply(cells_per_pert, function(n_cells) {
    simulate_one_pert(cells, guides = guides_per_pert, cells_per_pert = n_cells)
  })
  
  # split output into cre_perts and grna_perts
  cre_perts <- lapply(pert_status, FUN = "[[", 1)
  grna_perts <- lapply(pert_status, FUN = "[[", 2)
  
  # combine into 2 matrices
  cre_perts <- do.call(rbind, cre_perts)
  grna_perts <- do.call(rbind, grna_perts)
  
  # set rownames with enhancer and grna ids
  enh_ids <- paste0("enh", seq(length(num_cells_per_pert)))
  grna_ids <- paste0(rep(enh_ids, each = guides_per_pert), paste0("_g", seq_len(guides_per_pert)))
  rownames(cre_perts) <- enh_ids
  rownames(grna_perts) <- grna_ids
  
  # add to sce as perturbation status alt experiments
  altExp(sce, e = "cre_perts") <- SummarizedExperiment(assays = list(perts = cre_perts))
  altExp(sce, e = "grna_perts") <- SummarizedExperiment(assays = list(perts = grna_perts))
  
  return(sce)
  
}

# randomly draw perturbed cells and return a CRE and gRNA perturbation status matrix
simulate_one_pert <- function(cells, guides, cells_per_pert) {
  
  # draw number of cells per guide from a poisson with lambda = cells_per_guide
  n_cells <- rpois(guides, lambda = round(cells_per_pert / guides))
  
  # randomly draw cells perturbed by these guides
  pert_cells <- sample(seq_along(cells), size = sum(n_cells))

  # create CRE perturbations for these cells
  cre_perts <- matrix(rep(0, times = length(cells)), nrow = 1, dimnames = list(NULL, cells))
  cre_perts[, pert_cells] <- 1
  
  # create gRNA perturbations for these cells
  grna_perts <- matrix(0, nrow = guides, ncol = length(cells), dimnames = list(NULL, cells))
  pert_cells_per_grna <- split(pert_cells, f = rep(seq_len(guides), times = n_cells))
  for (grna in names(pert_cells_per_grna)) {
    grna_perts[as.integer(grna), pert_cells_per_grna[[grna]]] <- 1
  }
  
  # return perturbation status matrices
  return(list(cre_perts = cre_perts, grna_perts = grna_perts))
  
}


simulate_sample_perturbation_response <- function(pert, sce, effect_size, guide_sd, grna_target_data_frame) {
  # Initialize the pert object with the given pert and sce
  pert_object <- pert_input(pert, sce = sce, pert_level = "cre_perts")
  
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
  es_mat <- create_effect_size_matrix(grna_pert_status, pert_guides = pert_guides,
                                      gene_effect_sizes = effect_sizes, guide_sd = guide_sd)
  
  # Center effect sizes on specified gene-level effect sizes
  es_mat <- center_effect_size_matrix(es_mat, pert_status = pert_status, gene_effect_sizes = effect_sizes)
  es_mat_use = es_mat[, colnames(assay(pert_object, "counts"))]
  
  # Simulate Counts
  message("Simulating Counts")
  sim_counts <- sim_tapseq_sce(pert_object, effect_size_mat = es_mat_use)
  
  # Return the simulated response matrix as a sparse matrix
  simulated_response_matrix <- as(assay(sim_counts, "counts"), "RsparseMatrix")
  return(simulated_response_matrix)
}
