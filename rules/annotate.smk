# Dynamically resolve the correct bacon or meta-analysis input file
def get_file(wildcards):
    prefix = config["out_prefix"]
    assoc = config["association_variable"]
    out_dir = config["out_directory"]
    out_type = config["out_type"]

    if config["stratified_ewas"] == "yes":
        return f"{out_dir}{prefix}_{assoc}_ewas_meta_analysis_results_1.txt"
    else:
        return f"{out_dir}{prefix}_{assoc}_ewas_bacon_results{out_type}"

rule add_annotation:
    input: 
        annotation = config["annotation_manifest"],
        snp_key = config["snp_annotation"],
        in_file = get_file,
        script = "scripts/annotation.R"
    params:
        o_dir = config["out_directory"],
        strat = config["stratified_ewas"],
        assoc = config["association_variable"],
        o_type = config["out_type"],
        prefix = config["out_prefix"]
    output: 
        lambda wildcards: f"{config['out_directory']}{config['out_prefix']}_{config['association_variable']}_ewas_annotated_results{config['out_type']}"
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        Rscript {input.script} \
        --input-file {input.in_file} \
        --out-dir {params.o_dir} \
        --stratified {params.strat} \
        --assoc {params.assoc} \
        --out-type {params.o_type} \
        --annotation {input.annotation} \
        --snp-key {input.snp_key}
        """
