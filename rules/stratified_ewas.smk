rule stratify_data:
    input:
        script = "scripts/stratify.R",
        pheno_file = PHENO,
        methyl_file = MVALS
    params:
        strat_vars = ' '.join(STRAT_VARS),
        o_dir = OUT_DIR,
        n_threads = N_WORKERS
    output: 
        temp(expand(OUT_DIR + "{group}/{group}_pheno.fst", group = GROUPS)),
        temp(expand(OUT_DIR + "{group}/{group}_mvals.fst", group = GROUPS))
    conda:
        "../envs/ewas.yaml"    
    shell:
        """
        Rscript {input.script} \
        --pheno {input.pheno_file} \
        --methyl {input.methyl_file} \
        --stratify {params.strat_vars} \
        --out-dir {params.o_dir} \
        --threads {params.n_threads}
        """

rule run_ewas_stratified:
    input:
        script = "scripts/ewas.R",
        pheno_file = lambda wc: f"{OUT_DIR}{wc.group}/{wc.group}_pheno.fst",
        methyl_file = lambda wc: f"{OUT_DIR}{wc.group}/{wc.group}_mvals.fst"
    output:
        result = lambda wc: f"{OUT_DIR}{wc.group}/{wc.group}_{ASSOC}_ewas_results{OUT_TYPE}"
    log:
        lambda wc: f"log/{wc.group}_ewas.log"
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        export R_PROGRESSR_ENABLE=TRUE 
        Rscript {input.script} \
        --pheno {input.pheno_file} \
        --methyl {input.methyl_file} \
        --assoc {ASSOC} \
        --stratified {STRATIFIED} \
        --chunk-size {CHUNK_SIZE} \
        --processing-type {PROCESSING_TYPE} \
        --workers {N_WORKERS} \
        --out-dir {OUT_DIR}{wildcards.group}/ \
        --out-type {OUT_TYPE} \
        --out-prefix {wildcards.group} \
        > {log.logfile} 2>&1
        """

rule run_bacon_stratified:
    input:
        in_file = lambda wc: f"{OUT_DIR}{wc.group}/{wc.group}_{ASSOC}_ewas_results{OUT_TYPE}",
        script = "scripts/run_bacon.R"
    output:
        result = lambda wc: f"{OUT_DIR}{wc.group}/{wc.group}_{ASSOC}_ewas_bacon_results{OUT_TYPE}",
        trace = lambda wc: f"{OUT_DIR}{wc.group}/bacon_plots/{wc.group}_{ASSOC}_traces.jpg",
        post = lambda wc: f"{OUT_DIR}{wc.group}/bacon_plots/{wc.group}_{ASSOC}_posteriors.jpg",
        fit = lambda wc: f"{OUT_DIR}{wc.group}/bacon_plots/{wc.group}_{ASSOC}_fit.jpg",
        qq = lambda wc: f"{OUT_DIR}{wc.group}/bacon_plots/{wc.group}_{ASSOC}_qqs.jpg"
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        Rscript {input.script} \
        --input-file {input.in_file} \
        --out-dir {OUT_DIR}{wildcards.group}/ \
        --out-prefix {wildcards.group} \
        --out-type {OUT_TYPE}
        """

rule make_metal_script:
    input:
        script = "scripts/metal_cmd.sh",
        in_files = expand(OUT_DIR + "{group}/{group}_" + ASSOC + "_ewas_bacon_results" + OUT_TYPE, group=GROUPS)
    params:
        out_prefix = OUT_DIR + ASSOC + "_ewas_meta_analysis_results_"
    output:
        "scripts/meta_analysis_script.sh"
    shell:
        "sh {input.script} {input.in_files} {params.out_prefix}"

rule run_metal:
    input: 
        script = "scripts/meta_analysis_script.sh"
    output:
        meta_analysis_results
    singularity:
        "library://krferrier/metal/meta_analysis:metal"
    shell: 
        "metal {input.script}"
