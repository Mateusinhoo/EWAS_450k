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
