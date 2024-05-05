
# These rules are used to simulate the guide assignments given parameters in the config file
# Along with the simulation, dispersion estimates are calculated on the real data
    
# Create the simulated guide and cre assignments
rule simulate_guide_assignments:
  input:
    raw_counts = "resources/{sample}/raw_counts.rds"
  output:
    simulated_sce = "results/{sample}/simulated_sce.rds"
  params:
    num_cells_per_pert = config["simulate_guide_assignments"]["num_cells_per_pert"],
    num_guides_per_pert = config["simulate_guide_assignments"]["num_guides_per_pert"]
  log: 
    "results/{sample}/logs/simulate_guide_assignment.log"
  conda:
    "../envs/sceptre_power_simulations.yml"
  resources:
    mem = "64G",
    time = "2:00:00"
  script:
    "../scripts/simulate_guide_assignments.R"
    
# Fit dispersion estimates and calculate size_factors
rule fit_dispersions:
  input:
    raw_counts = "resources/{sample}/raw_counts.rds",
    simulated_sce = "results/{sample}/simulated_sce.rds"
  output:
    simulated_sce_disp = "results/{sample}/simulated_sce_disp.rds"
  params:
    size_factors = config["fit_dispersions"]["size_factors"],
    fit_type = config["fit_dispersions"]["fit_type"],
    disp_type = config["fit_dispersions"]["disp_type"]
  log: 
    "results/{sample}/logs/fit_dispersions.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "4:00:00"
  script:
    "../scripts/fit_negbinom_distr.R"
