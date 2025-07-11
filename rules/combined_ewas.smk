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
        o_prefix = OUT_PREFIX,
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
        --out-prefix {params.o_prefix} \
        """

rule bacon_correction:
    input:
        raw_results
    output:
        bacon_results,
        *bacon_plots  # these are the images created by ggsave()
    params:
        out_prefix = OUT_PREFIX,
        out_dir = OUT_DIR,
        out_type = OUT_TYPE
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
