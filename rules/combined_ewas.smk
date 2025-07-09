rule run_combined_ewas:
    input:
        script = "scripts/ewas.R",
        pheno_file = PHENO,
        methyl_file = MVALS
    params:
        assoc_var = ASSOC,
        stratified = STRATIFIED,
        cs = CHUNK_SIZE,
        pt = PROCESSING_TYPE,
        n_workers = N_WORKERS,
        o_dir = OUT_DIR,
        o_type = OUT_TYPE,
        o_prefix = OUT_PREFIX  
    output: 
        raw_results
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

rule run_bacon:
    input:
        in_file = raw_results,
        script = "scripts/run_bacon.R"
    params:
        o_dir = OUT_DIR,
        o_type = OUT_TYPE,
        o_prefix = OUT_PREFIX + "_" + ASSOC
    output:
        f"{OUT_DIR}{OUT_PREFIX}_{ASSOC}_ewas_bacon_results{OUT_TYPE}",
        f"{OUT_DIR}bacon_plots/{OUT_PREFIX}_{ASSOC}_traces.jpg",
        f"{OUT_DIR}bacon_plots/{OUT_PREFIX}_{ASSOC}_posteriors.jpg",
        f"{OUT_DIR}bacon_plots/{OUT_PREFIX}_{ASSOC}_fit.jpg",
        f"{OUT_DIR}bacon_plots/{OUT_PREFIX}_{ASSOC}_qqs.jpg"
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        Rscript {input.script} \
        --input-file {input.in_file} \
        --out-dir {params.o_dir} \
        --out-prefix {params.o_prefix} \
        --out-type {params.o_type}
        """
