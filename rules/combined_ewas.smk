rule run_combined_ewas:
    input:
        script = "scripts/ewas.R",
        pheno_file = config["pheno"],
        methyl_file = config["mvals"]
    params:
        assoc_var = config["association_variable"],
        stratified = config["stratified_ewas"],
        cs = config["chunk_size"],
        pt = config["processing_type"],
        n_workers = config["workers"],
        o_dir = config["out_directory"],
        o_type = config["out_type"],
        o_prefix = config["out_prefix"]
    output: 
        config["out_directory"] + config["out_prefix"] + "_" + config["association_variable"] + "_ewas_results" + config["out_type"]
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        Rscript -e "Sys.setenv(R_PROGRESSR_ENABLE='TRUE'); source('{input.script}')" \
        --pheno {input.pheno_file} \
        --methyl {input.methyl_file} \
        --assoc {params.assoc_var} \
        --stratified {params.stratified} \
        --chunk-size {params.cs} \
        --processing-type {params.pt} \
        --workers {params.n_workers} \
        --out-dir {params.o_dir} \
        --out-type {params.o_type} \
        --out-prefix {params.o_prefix}
        """

rule bacon_correction:
    input:
        config["out_directory"] + config["out_prefix"] + "_" + config["association_variable"] + "_ewas_results" + config["out_type"]
    output:
        config["out_directory"] + config["out_prefix"] + "_" + config["association_variable"] + "_ewas_bacon_results" + config["out_type"],
        directory(config["out_directory"] + "/bacon_plots")
    params:
        out_prefix = config["out_prefix"],
        out_dir = config["out_directory"],
        out_type = config["out_type"]
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        mkdir -p {params.out_dir}/bacon_plots
        Rscript scripts/run_bacon.R \
            --input-file {input} \
            --out-dir {params.out_dir} \
            --out-prefix {params.out_prefix} \
            --out-type {params.out_type}
        """
