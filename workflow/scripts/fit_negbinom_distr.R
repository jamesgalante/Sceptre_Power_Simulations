## Fit negative binomial distributions for each gene using DESeq2

# Saving image for debugging
save.image("RDA_objects/fit_dispersions.rda")
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# opening log file to collect all messages, warnings and errors
message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages and functions
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(DESeq2)
  source(file.path(snakemake@scriptdir, "R_functions/power_simulations_fun.R"))
})

# load prepared input data stored in SingleCellExperiment object
message("Loading input data.")
response_matrix <- readRDS(snakemake@input$raw_counts)
colnames(response_matrix) <- NULL
simulated_sce <- readRDS(snakemake@input$simulated_sce)

# Calculate total_umis and detected_genes for Deseq2 object creation
coldata <- data.frame(
  total_umis = colSums(response_matrix),
  detected_genes = colSums(response_matrix > 0)
)

# fit negative binomial distributions to estimate gene-level dispersion
message("Estimate dispersion using DESeq2:")
sce <- fit_negbinom_deseq2(response_matrix,
                           simulated_sce,
                           coldata,
                           size_factors = "poscounts",
                           fit_type = "parametric",
                           disp_type = "dispersion")


# save output sce to file
saveRDS(sce, file = snakemake@output$simulated_sce_disp)

# close log file connection
sink()
sink(type = "message")
close(log)





