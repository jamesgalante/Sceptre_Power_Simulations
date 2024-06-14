# Script to simulate guide assignments given parameters in config file

### CREATE DEBUG FILES =======================================================

# Saving image for debugging
save.image("RDA_objects/simulate_guide_assignments.rda")
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
  source(file.path(snakemake@scriptdir, "R_functions/simulate_perturbations.R"))
  source(file.path(snakemake@scriptdir, "R_functions/power_simulations_fun.R"))
})

# Download Sceptre
message("Installing Sceptre")
library(devtools)
devtools::install_github("katsevich-lab/sceptre")
message("Sceptre Installation Complete")
library(sceptre)

# Load in input and the response matrix
message("Loading input and response matrix")
response_matrix <- readRDS(snakemake@input$raw_counts)

# Load in snakemake parameters
message("Loading parameters")
num_cells_per_pert <- snakemake@params$num_cells_per_pert
num_guides_per_pert <- snakemake@params$num_guides_per_pert
num_cells <- snakemake@params$num_cells

### SIMULATE GUIDE COUNTS ====================================================

message("Creating the sce object")
# Estimate how many cells we need in our sce object
if (is.null(num_cells)) {
  num_cells <- num_cells_per_pert[length(num_cells_per_pert)] * 5 + 6000 # The 6000 is to make sure there are enough n_ctrl cells but this calculation is arbitrary otherwise
  message(paste("Estimated number of cells:", num_cells))
}

# Create the raw matrix with the new number of cells
counts <- matrix(0, nrow=nrow(response_matrix), ncol=num_cells)
colnames(counts) <- paste("cell", 1:num_cells, sep="")
rownames(counts) <- rownames(response_matrix)

# Create the sce object
sce <- SingleCellExperiment(assays=list(counts=counts))

# Simulate perturbations
message("Simulating perturbations")
sce <- simulate_perturbations(sce, 
                              cells_per_pert=num_cells_per_pert, 
                              guides_per_pert=num_guides_per_pert)


# Fill in row data names for grna_perts
rowData(altExp(sce, "grna_perts"))$name <- rownames(altExp(sce, "grna_perts"))
rowData(altExp(sce, "grna_perts"))$target_name <- sub("_.*", "", rownames(altExp(sce, "grna_perts")))

# save simulation output
message("Saving output to file.")
saveRDS(sce, file = snakemake@output$simulated_sce)

# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)