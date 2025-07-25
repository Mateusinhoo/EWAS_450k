default_container: "docker://rocker/tidyverse:4.3.1"

import pandas as pd
from helper_fxns import generate_observed_combinations
configfile: "config.yml"

#----SET VARIABLES----#
## EWAS VARIABLES
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

# DMR VARIABLES
DMR = config["dmr_analysis"]
ANNO = config["genome_build"]

if DMR == "yes":
    MIN_P = config["min_pvalue"]
    WIN_SZ = config["window_size"]
    REGION_FILTER = config["region_filter"]

# Stratified EWAS
if STRATIFIED == "yes":
    GROUPS = generate_observed_combinations(
        df=pd.read_csv(config["pheno"]),
        stratify_cols=config["stratify_variables"]
    )
else:
    GROUPS = []

#---- INPUT & OUTPUT FILES ----#
# Final output results, stratified or not
annotated_results = OUT_DIR + OUT_PREFIX + "_" + ASSOC + "_ewas_annotated_results" + OUT_TYPE
manhattan_qq_plot = OUT_DIR + OUT_PREFIX + "_" + ASSOC + "_ewas_manhattan_qq_plots.jpg"

# Combined (not stratified) EWAS outputs
raw_results = OUT_DIR + OUT_PREFIX + "_" + ASSOC + "_ewas_results" + OUT_TYPE
bacon_results = OUT_DIR + OUT_PREFIX + "_" + ASSOC + "_ewas_bacon_results" + OUT_TYPE
bacon_plots = expand(OUT_DIR + "bacon_plots/" + OUT_PREFIX + "_" + ASSOC + "_{plot}.jpg", plot=PLOTS)

# Stratified EWAS outputs
strat_raw_results = expand(OUT_DIR + "{group}/{group}_" + OUT_PREFIX + "_" + ASSOC + "_ewas_results" + OUT_TYPE, group=GROUPS)
strat_bacon_results = expand(OUT_DIR + "{group}/{group}_" + OUT_PREFIX + "_" + ASSOC + "_ewas_bacon_results" + OUT_TYPE, group=GROUPS)
strat_bacon_plots = expand(OUT_DIR + "{group}/bacon_plots/{group}_" + OUT_PREFIX + "_" + ASSOC + "_{plot}.jpg", group=GROUPS, plot=PLOTS)
meta_analysis_results = OUT_DIR + OUT_PREFIX + "_" + ASSOC + "_ewas_meta_analysis_results_1.txt"

# DMR outputs
results_bed = OUT_DIR + OUT_PREFIX + "_" + ASSOC + "_ewas_annotated_results.bed"
dmr_acf = OUT_DIR + "dmr/" + OUT_PREFIX + "_" + ASSOC + "_ewas.acf.txt"
dmr_args = OUT_DIR + "dmr/" + OUT_PREFIX + "_" + ASSOC + "_ewas.args.txt"
dmr_fdr = OUT_DIR + "dmr/" + OUT_PREFIX + "_" + ASSOC + "_ewas.fdr.bed.gz"
dmr_regions = OUT_DIR + "dmr/" + OUT_PREFIX + "_" + ASSOC + "_ewas.regions.bed.gz"
dmr_slk = OUT_DIR + "dmr/" + OUT_PREFIX + "_" + ASSOC + "_ewas.slk.bed.gz"

dmr_infile = [results_bed]
dmr_outfiles = [dmr_acf, dmr_args, dmr_fdr, dmr_regions, dmr_slk]

#----RULE ALL----#
rule all:
    input:
        "results/.placeholder",
        PHENO,
        MVALS,
        # EWAS results
        OUT_DIR + OUT_PREFIX + "_" + ASSOC + "_ewas_results" + OUT_TYPE,
        OUT_DIR + OUT_PREFIX + "_" + ASSOC + "_ewas_bacon_results" + OUT_TYPE,
        expand(OUT_DIR + "bacon_plots/" + OUT_PREFIX + "_" + ASSOC + "_{plot}.jpg", plot=PLOTS),
        annotated_results,
        manhattan_qq_plot

#----INCLUDE RULES----#
include: "rules/combined_ewas.smk"
include: "rules/annotate.smk"
include: "rules/plots.smk"
include: "rules/setup.smk"
