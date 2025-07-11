rule make_results_dir:
    output:
        "results/.placeholder"
    shell:
        "mkdir -p results && touch {output}"
