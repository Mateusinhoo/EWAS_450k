rule make_results_dir:
    output:
        directory(config["out_directory"])
    shell:
        "mkdir -p {output}"

rule preprocess_idat:
    input:
        idat_dir = "data/idat/",          # Folder with .idat files
        sample_sheet = "data/sample_sheet.csv"
    output:
        mvals = "data/mvals_minfi.csv"
    script:
        "scripts/preprocess_idat_minfi.R"
