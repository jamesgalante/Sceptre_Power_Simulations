# Sceptre-based Power Simulations

## About The Project

The **Sceptre-based Power Simulations** pipeline is designed to help determine the number of cells needed in CRISPRi or CRISPRko experiments to achieve a desired statistical power (e.g., 80%) for detecting a specific percentage decrease in gene expression. By simulating different effect sizes and cell numbers, this tool provides valuable insights for experimental design and helps determine the feasibility of detecting small changes in gene expression.

This pipeline leverages the **sceptre** R package, which offers a statistically rigorous and scalable approach for single-cell CRISPR screen data analysis. It calculates mean expression and dispersion from your input data, simulates a specified percentage decrease in gene expression, and computes the statistical power to detect that change across varying numbers of cells.

## Installation

To get started with the Sceptre-based Power Simulations pipeline, follow these steps:

1. **Clone the Repository**

   ```bash
   git clone https://github.com/jamesgalante/Sceptre_Power_Simulations.git
   cd Sceptre_Power_Simulations
   ```

2. **Create a Conda Environment with Snakemake**

   It's recommended to create a new conda environment specifically for this project. We'll install Snakemake version `7.32.4`.

   ```bash
   conda create -n sceptre_power_sim python=3.9
   conda activate sceptre_power_sim
   conda install -c bioconda snakemake=7.32.4
   ```

   *Note:* All other required packages will be dynamically downloaded when the Snakemake pipeline is run.

## Usage

### Input Data Preparation

The pipeline requires a raw counts matrix in the form of a `.rds` file (created with saveRDS). This file should be a `dgRMatrix` format. An example of the expected format can be found in `resources/test_data/raw_counts.rds`.

If you wish to visualize the power simulations results by TPM instead of UMIs, provide a TPM file formatted as in `resources/tpm_per_gene.tsv`. The TPM file is only used to match your input genes with the power simulation results during plotting. If you are designing a screen based on TPM values, this might be a good option. Since single-cell UMIs may differ from TPM calculated values in bulk screens, the resulting plots may look more jagged.

### Configuring the Pipeline

Before running the pipeline, you need to set up the configuration file `config.yaml`. This file controls various parameters of the simulation and analysis.

Here is an example of the `config.yaml` file:

```yaml
# Main Config File Format for Sceptre-based Power Simulations Pipeline

samples:
  your_sample_name

simulate_guide_assignments:
  num_cells_per_pert: [50, 75, 100, 150, 250, 500, 750, 1000, 1400, 1800, 2500, 4000, 7500, 10000]

sceptre_power_analysis:
  effect_size: 0.15
  reps: 20

compute_power_from_simulations:
  pval_adj_thresh: 0.1
  positive_proportion: 0.05

visualize_power_results:
  tpm_per_gene: "resources/tpm_per_gene.tsv"  # Set to NULL if not using TPM
  gene_format: "ensembl"  # Must be "ensembl" or "symbol"
```

#### Configuration Parameters Explained

- **samples**: A list containing the names of the samples you wish to analyze. Each sample name should correspond to a directory in `resources/` that contains the `raw_counts.rds` file. Replace `your_sample_name` with your actual sample name.

- **simulate_guide_assignments**:
  - **num_cells_per_pert**: A list of numbers representing the different counts of cells per perturbation you want to simulate. The pipeline will assess the statistical power for each specified cell count.

- **sceptre_power_analysis**:
  - **effect_size**: The percent decrease in gene expression you wish to simulate. For example, `0.15` represents a 15% decrease.
  - **reps**: The number of simulation replicates to run for each condition. More replicates reduce noise but increase computational time. A value of `20` is typically sufficient.

- **compute_power_from_simulations**:
  - **pval_adj_thresh**: The adjusted p-value threshold for determining statistical significance after applying Benjamini-Hochberg (BH) correction. The default Sceptre recommendation is `0.1`.
  - **positive_proportion**: The expected proportion of true positives (genes with actual expression changes) in your data. This parameter is important because the BH correction assumes a certain proportion of null hypotheses. This number is hard to estimate before running an experiment, but for CRISPR screens, the proportion of positive tests is typically low. A value of `0.05` corresponds to 5% of the resulting tests being positives.

- **visualize_power_results**:
  - **tpm_per_gene**: Path to the TPM per gene file. If you wish to use UMI counts instead, set this to `NULL`. See the example `resources/tpm_per_gene.tsv` for formatting requirements.
  - **gene_format**: Specifies the format of gene names used in your data. Must be either `"ensembl"` or `"symbol"` or refer to a specific column in the TPM file.

### Running the Pipeline

To execute the pipeline, run the following commands from the root directory of the project:

1. **Dry Run**

   It's good practice to perform a dry run to ensure that the pipeline is configured correctly and all files are in place.

   ```bash
   snakemake --use-conda all -np
   ```

   This command will show you what steps the pipeline will perform without actually executing them.

2. **Full Run**

   If the dry run looks correct, execute the pipeline:

   ```bash
   snakemake --use-conda all
   ```

   *Note:* The `--use-conda` flag tells Snakemake to create and use the necessary conda environments for each rule.
   *Note:* See the relevant documentation [Snakemake Documentation](https://snakemake.readthedocs.io/en/stable/index.html) to understand what flags might be necessary. If using slurm, it might be easiest to create a snakemake profile in `~/.config/snakemake/profile_name/config.yml`. Here is an example profile:
   
   ```yaml
   jobs: 100 
   slurm: True 
   retries: 1
   use-conda: True 
   notemp: True 
   default-resources: 
       - slurm_account=slurm_account_name
       - slurm_partition=slurm_account_name,owners,normal
       - runtime="6h"
       - slurm_extra="--nice"
   ```
   The pipeline can then be run with `snakemake --profile profile_name all`

#### Troubleshooting

- If you encounter issues, you can test the pipeline with the provided test data. Rename or delete the `results/test_data/` directory and run the pipeline again to ensure everything works as expected.

- Ensure that your `config.yaml` file is properly configured and that the `raw_counts.rds` object matches the format of the `test_data`.


## How It Works

The pipeline simulates gene expression data to assess the statistical power of detecting specified effect sizes across varying numbers of cells. It begins by calculating mean expression and dispersion for each gene from your input data. Cells are then assigned to control or perturbation groups, simulating a decrease in expression based on the specified effect size. Using the sceptre package, differential expression analysis is performed, and p-values are adjusted using the Benjamini-Hochberg procedure, considering the expected proportion of true positives. Finally, the pipeline computes the power for each condition and generates visualizations to help interpret the results.

## Notes

- **Data Format**: Ensure that your input data is correctly formatted. The gene names should match between your raw counts matrix and the TPM file (if used).

- **Sceptre Package**: The sceptre package is maintained by another team. This pipeline utilizes it for its efficiency and statistical robustness in single-cell CRISPR screen analysis. For more details, refer to their documentation.

- **Understanding FDR Correction**: The Benjamini-Hochberg procedure assumes a certain proportion of true positives. The `positive_proportion` parameter allows you to adjust this assumption based on your experimental context.

- **TO DO**: Use dispersion estimates calculated by Sceptre

## References

- **Sceptre Package Documentation**: [Hands-On Single-Cell CRISPR Screen Analysis](https://timothy-barry.github.io/sceptre-book/)

---
