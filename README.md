# Sceptre_Power_Simulations
Statistical power simulator with customizable inputs


To run this pipeline all you need is the response_matrix from a sceptre object, i.e. a raw counts matrix
stored in an rds file as a dgRMatrix, which can be seen as follows:

"""
> str(response_matrix)
Formal class 'dgRMatrix' [package "Matrix"] with 6 slots
  ..@ p       : int [1:527] 0 22616 29856 30002 38424 42120 46117 46358 50797 52091 ...
  ..@ j       : int [1:3986055] 0 1 2 3 4 5 6 7 8 9 ...
  ..@ Dim     : int [1:2] 526 24425
  ..@ Dimnames:List of 2
  .. ..$ : chr [1:526] "ENSG00000069275" "ENSG00000117222" "ENSG00000117266" "ENSG00000117280" ...
  .. ..$ : NULL
  ..@ x       : num [1:3986055] 4 6 10 9 4 13 8 8 10 8 ...
  ..@ factors : list()
"""

The example data is from this Sceptre script, which can be found in the Sceptre documentation.
For downloading Sceptre please also see Sceptre documentation: "https://timothy-barry.github.io/sceptre-book/"

# Load sceptre and sceptredata
library(sceptre)
library(sceptredata)
# load the data, creating a sceptre_object
directories <- paste0(
  system.file("extdata", package = "sceptredata"),
  "/highmoi_example/gem_group_", 1
)
sceptre_object <- import_data_from_cellranger(
  directories = directories,
  moi = "high",
  grna_target_data_frame = data(grna_target_data_frame_highmoi)
)
file_path <- "resources/test_data/raw_counts.rds"
saveRDS(get_response_matrix(sceptre_object), file_path)
  
  
To run everything:
  Run `snakemake --profile slurm_larsms all -np` in a tmux session where "slurm_larsms" is a snakemake profile stored
  in "~/.config/snakemake/slurm_larsms/config.yaml"
  