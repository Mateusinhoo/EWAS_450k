default_container: "docker://rocker/tidyverse:4.3.1"

import yaml
import pandas as pd
from helper_fxns import generate_observed_combinations

with open("config_all.yml", "r") as f:
    configs = yaml.safe_load(f)["configs"]

SAMPLES = [c.replace(".yml", "").replace("config_", "") for c in configs]

# Define plots
PLOTS = ["traces", "posteriors", "fit", "qqs"]

# Build all expected output paths dynamically
def get_outputs(sample):
    with open(f"config_{sample}.yml", "r") as f:
        cfg = yaml.safe_load(f)

    out_prefix = cfg["out_prefix"]
    out_dir = cfg["out_directory"]
    assoc = cfg["association_variable"]
    out_type = cfg["out_type"]
    stratified = cfg["stratified_ewas"]
    dmr = cfg.get("dmr_analysis", "no")

    # Basic files
    files = [
        cfg["pheno"],
        cfg["mvals"],
        f"{out_dir}{out_prefix}_{assoc}_ewas_results{out_type}",
        f"{out_dir}{out_prefix}_{assoc}_ewas_bacon_results{out_type}",
        f"{out_dir}{out_prefix}_{assoc}_ewas_annotated_results{out_type}",
        f"{out_dir}{out_prefix}_{assoc}_ewas_manhattan_qq_plots.jpg",
    ]

    files += [f"{out_dir}bacon_plots/{out_prefix}_{assoc}_{plot}.jpg" for plot in PLOTS]

    if stratified == "yes":
        df = pd.read_csv(cfg["pheno"])
        groups = generate_observed_combinations(df, cfg["stratify_variables"])
        for g in groups:
            files += [
                f"{out_dir}{g}/{g}_{out_prefix}_{assoc}_ewas_results{out_type}",
                f"{out_dir}{g}/{g}_{out_prefix}_{assoc}_ewas_bacon_results{out_type}",
            ]
            files += [f"{out_dir}{g}/bacon_plots/{g}_{out_prefix}_{assoc}_{plot}.jpg" for plot in PLOTS]
        files.append(f"{out_dir}{out_prefix}_{assoc}_ewas_meta_analysis_results_1.txt")

    if dmr == "yes":
        files += [
            f"{out_dir}{out_prefix}_{assoc}_ewas_annotated_results.bed",
            f"{out_dir}dmr/{out_prefix}_{assoc}_ewas.acf.txt",
            f"{out_dir}dmr/{out_prefix}_{assoc}_ewas.args.txt",
            f"{out_dir}dmr/{out_prefix}_{assoc}_ewas.fdr.bed.gz",
            f"{out_dir}dmr/{out_prefix}_{assoc}_ewas.regions.bed.gz",
            f"{out_dir}dmr/{out_prefix}_{assoc}_ewas.slk.bed.gz",
        ]

    return files

# Rule all combines all outputs from all configs
rule all:
    input:
        expand(get_outputs, sample=SAMPLES)

include: "rules/combined_ewas.smk"
include: "rules/stratified_ewas.smk"
include: "rules/annotate.smk"
include: "rules/plots.smk"
# include: "rules/dmr.smk"
