# These rule split the target-response pairs file into batches to be processed in parallel
# These outputs are then recombined after power simulations

# Create and split the discovery pairs file
rule split_target_response_pairs:
  input:
    simulated_sce_disp = "results/{sample}/simulated_sce_disp.rds"
  output:
    discovery_pairs = "resources/{sample}/discovery_pairs.txt",
    splits = expand("resources/{{sample}}/pair_splits/discovery_pairs_split_{split}.txt",split=range(1,config['split_target_response_pairs']['batches'] + 1)),
    grna_target_data_frame = "resources/{sample}/grna_target_data_frame.txt"
  params:
    batches = config["split_target_response_pairs"]['batches']
  log: "results/{sample}/logs/split_target_response_pairs.log"
  conda:
    "../envs/sceptre_power_simulations.yml"
  resources:
    mem = "16G",
    time = "1:00:00"
  script:
    "../scripts/split_target_response_pairs.R"
    
    
# Let's simulate counts and create a proxy sceptre object that can be passed into each power simulation
# We do this, so we only have to change the covariates when the response_matrix is simulated in each rep
rule create_simulated_sceptre_object:
  input:
    simulated_sce_disp = "results/{sample}/simulated_sce_disp.rds",
    grna_target_data_frame = "resources/{sample}/grna_target_data_frame.txt",
    discovery_pairs = "resources/{sample}/discovery_pairs.txt"
  output:
    simulated_sceptre_object = "resources/{sample}/simulated_sceptre_object.rds"
  log: "results/{sample}/logs/create_simulated_sceptre_object.log"
  conda:
    "../envs/sceptre_power_simulations.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/create_simulated_sceptre_object.R"


# Run the power simulation with sceptre for each split
rule sceptre_power_analysis:
  input:
    # Maybe change these names to match the output names of the rules that lead up to it
    simulated_sceptre_object	= "resources/{sample}/simulated_sceptre_object.rds",
    simulated_sce_disp	= "results/{sample}/simulated_sce_disp.rds",
    discovery_pairs_split = "resources/{sample}/pair_splits/discovery_pairs_split_{split}.txt",
    grna_target_data_frame = "resources/{sample}/grna_target_data_frame.txt",
    raw_counts = "resources/{sample}/raw_counts.rds"
  output:
    "results/{sample}/power_analysis_split/power_analysis_output_{split}.tsv"
  params:
    effect_size = config["sceptre_power_analysis"]["effect_size"],
    reps = config["sceptre_power_analysis"]["reps"],
    guide_sd = config["sceptre_power_analysis"]["guide_sd"]
  log: "results/{sample}/logs/sceptre_power_analysis_{split}.log"
  conda:
    "../envs/sceptre_power_simulations.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/sceptre_power_analysis.R"


# Combine the split outputs of the power analysis
rule combine_sceptre_power_analysis:
 input:
   expand("results/{{sample}}/power_analysis_split/power_analysis_output_{split}.tsv", split=range(1, config["split_target_response_pairs"]["batches"] + 1))
 output:
   "results/{sample}/power_analysis_output.tsv"
 log: "results/{sample}/logs/combine_sceptre_power_analysis.log"
 conda:
   "../envs/sceptre_power_simulations.yml"
 resources:
   mem = "32G",
   time = "2:00:00"
 script:
   "../scripts/combine_sceptre_power_analysis.R"
