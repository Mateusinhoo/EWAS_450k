rule make_results_dir:
    output:
        directory(config["out_directory"])
    shell:
        "mkdir -p {output}"

rule preprocess_idat:
    input:
        idat_dir = "/Volumes/LaCie/data",    #Raw Illumina data files
        sample_sheet = "data/sample_sheet.csv"
    output:
        mvals = "data/kidney_mvals.csv"
    script:
        "scripts/preprocess_idat_minfi.R"
