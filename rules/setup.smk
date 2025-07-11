rule make_results_dir:
    output:
        "results/"
    shell:
        "mkdir -p {output}"
