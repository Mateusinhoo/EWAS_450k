import yaml
import pandas as pd
from helper_fxns import generate_observed_combinations

configs = yaml.safe_load(open("config_all.yml"))["configs"]

rules_to_run = []

for config_path in configs:
    config = yaml.safe_load(open(config_path))
    
    # Pull variables from each config
    PHENO = config["pheno"]
    MVALS = config["mvals"]
    ASSOC = config["association_variable"]
    STRATIFIED = config["stratified_ewas"]
    STRAT_VARS = config["stratify_variables"]
    CHUNK_SIZE = config["chunk_size"]
    PROCESSING_TYPE = config["processing_type"]
    N_WORKERS = config["workers"]
    OUT_DIR = config["out_directory"]
    OUT_TYPE = config["out_type"]
    OUT_PREFIX = config["out_prefix"]
    ANNOTATION_MANIFEST = config["annotation_manifest"]
    SNP_ANNOTATION = config["snp_annotation"]
    PLOTS = ["traces", "posteriors", "fit", "qqs"]

    DMR = config.get("dmr_analysis", "no")

    if STRATIFIED == "yes":
        GROUPS = generate_observed_combinations(
            df=pd.read_csv(PHENO),
            stratify_cols=STRAT_VARS
        )
    else:
        GROUPS = []

    # Collect the input files per config
    rules_to_run += [
        PHENO,
        MVALS,
        *(expand(OUT_DIR + "{group}/{group}_" + OUT_PREFIX + "_" + ASSOC + "_ewas_results" + OUT_TYPE, group=GROUPS) if STRATIFIED == "yes" else [OUT_DIR + OUT_PREFIX + "_" + ASSOC + "_ewas_results" + OUT_TYPE]),
        *(expand(OUT_DIR + "{group}/{group}_" + OUT_PREFIX + "_" + ASSOC + "_ewas_bacon_results" + OUT_TYPE, group=GROUPS) if STRATIFIED == "yes" else [OUT_DIR + OUT_PREFIX + "_" + ASSOC + "_ewas_bacon_results" + OUT_TYPE]),
        *(expand(OUT_DIR + "{group}/bacon_plots/{group}_" + OUT_PREFIX + "_" + ASSOC + "_{plot}.jpg", group=GROUPS, plot=PLOTS) if STRATIFIED == "yes" else expand(OUT_DIR + "bacon_plots/" + OUT_PREFIX + "_" + ASSOC + "_{plot}.jpg", plot=PLOTS)),
        OUT_DIR + OUT_PREFIX + "_" + ASSOC + "_ewas_annotated_results" + OUT_TYPE,
        OUT_DIR + OUT_PREFIX + "_" + ASSOC + "_ewas_manhattan_qq_plots.jpg",
    ]

# Set the final target
rule all:
    input:
        rules_to_run

include: "rules/combined_ewas.smk"
include: "rules/stratified_ewas.smk"
include: "rules/annotate.smk"
include: "rules/plots.smk"
