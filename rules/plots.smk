rule plot_results:
    input: 
        in_file = config["out_directory"] + config["out_prefix"] + "_" + config["association_variable"] + "_ewas_annotated_results" + config["out_type"],
        script = "scripts/plots.R"
    params:
        o_dir = config["out_directory"],
        strat = config["stratified_ewas"],
        assoc = config["association_variable"],
        prefix = config["out_prefix"]
    output: 
       config["out_directory"] + config["out_prefix"] + "_" + config["association_variable"] + "_ewas_manhattan_qq_plots.jpg"
    conda:
        "../envs/ewas.yaml"
    shell:
        """
        Rscript {input.script} \
        --input-file {input.in_file} \
        --out-dir {params.o_dir} \
        --stratified {params.strat} \
        --assoc {params.assoc} \
        --out-prefix {params.prefix}
        """
