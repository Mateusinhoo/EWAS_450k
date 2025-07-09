# Manually construct expected input/output filenames
ewas_input_file = config["out_directory"] + config["out_prefix"] + "_" + config["association_variable"]
if config["stratified_ewas"] == "yes":
    ewas_input_file += "_ewas_meta_analysis_results_1.txt"
else:
    ewas_input_file += "_ewas_bacon_results" + config["out_type"]

annotated_output_file = config["out_directory"] + config["out_prefix"] + "_" + config["association_variable"] + "_ewas_annotated_results" + config["out_type"]

rule add_annotation:
    input: 
        annotation = config["annotation_manifest"],
        snp_key = config["snp_annotation"],
        in_file = ewas_input_file,
        script = "scripts/annotation.R"
    params:
        o_dir = config["out_directory"],
        strat = config["stratified_ewas"],
        assoc = config["association_variable"],
        o_type = config["out_type"]
    output: 
        annotated_output_file
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
