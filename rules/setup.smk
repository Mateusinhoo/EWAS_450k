rule make_results_dir:
    output:
        directory(OUT_DIR)
    shell:
        "mkdir -p {output}"
