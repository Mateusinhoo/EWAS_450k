def get_file(wildcards):
    if config["stratified_ewas"] == "yes":
        in_file = meta_analysis_results
    else:
        in_file = bacon_results
    return(in_file)

rule add_annotation:
    input: 
        annotation = config["annotation_manifest"],
        snp_key = config["snp_annotation"],
        in_file = get_file,
        script = "scripts/annotation.R"
    params:
        o_dir = OUT_DIR,
        strat = STRATIFIED,
        assoc = ASSOC,
        o_type = OUT_TYPE
    output: 
        annotated_results
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
