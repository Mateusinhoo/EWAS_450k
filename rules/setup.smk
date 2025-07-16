rule make_results_dir:
    output:
        directory(config["out_directory"])
    shell:
        "mkdir -p {output}"
