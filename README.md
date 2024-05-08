# Sceptre_Power_Simulations
Statistical power simulator with customizable inputs

To use
1. Create a `raw_counts.rds` file and store in `resources/{sample_name}/raw_counts.rds`. See example dataset for formatting
2. In the config file configure the range of cells per perturbation you want to simulate. This can be done with the `num_cells_per_pert` element. It's recommended to have the number of batches equal the length of the `num_cells_per_pert` list, so everything runs in parallel.
3. Run `snakemake {insert flags here} all -np`. If everything looks good, run without the `-np` flag
