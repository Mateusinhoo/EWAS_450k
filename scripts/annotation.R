# Script for running adding annotation data to ewas results
# Import libraries
suppressPackageStartupMessages({
    library(argparse)
    library(dplyr)
    library(data.table)
})
# Define command line arguments
parser <- argparse::ArgumentParser(description="Script for adding annotation data to ewas results")
parser$add_argument('--input-file', '-i',
                    required=TRUE,
                    help="Path to ewas results data")
parser$add_argument('--out-dir',
                    required=TRUE,
                    help="Path to output directory")
parser$add_argument('--stratified',
                    choices=c("yes", "no"), 
                    default="no",
                    help="Results from a stratified analysis: yes or no")
parser$add_argument("--assoc", 
                    required=TRUE,
                    type="character", 
                    nargs=1, 
                    help="Association variable EWAS was performed with.")
parser$add_argument('--out-type', 
                    type="character",
                    choices=c(".csv", ".csv.gz"), 
                    nargs="?",                    
                    const=".csv",
                    default=".csv",  
                    help="Output file type: CSV or CSV.GZ")

# parse arguments
args <- parser$parse_args()
results <- args$input_file
out_dir <- args$out_dir
stratified <- args$stratified
assoc <- args$assoc
out_type <- args$out_type

# Read in EWAS summary statistics
ewas <- fread(results)

# Ensure required Bioconductor annotation package is installed
# if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg38", quietly = TRUE)) {
#     if (!requireNamespace("BiocManager", quietly = TRUE))
#         install.packages("BiocManager")
#     BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg38", ask = FALSE)
# }

# Install bioconductor-qcewas
# if (!requireNamespace("QCEWAS", quietly = TRUE)) {
#     if (!requireNamespace("devtools", quietly = TRUE))
#         install.packages("devtools")
#     devtools::install_github("WaldronLab/QCEWAS", ask = FALSE)
# }

# Load annotation from TSV files
# ann_locations <- fread("annotation_files/EPIC_hg38.tsv.gz")
# ann_genes <- fread("annotation_files/EPIC_snp_key.tsv.gz")

# Make sure column names are set correctly
# colnames(ann_locations)[1] <- "cpgid"  # assuming first column is CpG ID
# colnames(ann_genes)[1] <- "cpgid"

# Merge them together
# annotation <- left_join(ann_locations, ann_genes, by = "cpgid")

# Convert Locations to a data frame with cpgid column
# ann_locations <- as.data.frame(Locations)
# ann_locations$cpgid <- rownames(ann_locations)

# Optionally add UCSC gene info (e.g., gene names)
# ann_genes <- as.data.frame(Other$UCSC_RefGene_Name)
# colnames(ann_genes) <- "gene"
# ann_genes$cpgid <- rownames(ann_genes)

# Merge into a final annotation table
# annotation <- dplyr::left_join(ann_locations, ann_genes, by = "cpgid")

# if(stratified=="no"){
#     ewas <- left_join(ewas, annotation, by = "cpgid")
#     ewas <- ewas[order(ewas$bacon.pval),]
# } else{
#         ewas <- left_join(ewas, annotation, by = c("MarkerName"="cpgid"))
#         ewas <- ewas %>% dplyr::select(-Allele1, -Allele2)
#         ewas$"P-value" <- as.numeric(ewas$"P-value")
#         ewas <- ewas[order(ewas$"P-value"),]
# }
# file_name <- paste0(out_dir, assoc, "_ewas_annotated_results", out_type)
# fwrite(ewas, file = file_name)

# Load local annotation files
ann_locations <- fread("annotation_files/EPIC_hg38.tsv.gz")  # CpG positions
ann_genes <- fread("annotation_files/EPIC_snp_key.tsv.gz")   # Gene annotations

# Ensure CpG ID columns match
colnames(ann_locations)[which(colnames(ann_locations) == "probeID")] <- "cpgid"
colnames(ann_genes)[which(colnames(ann_genes) == "probeID")] <- "cpgid"


# Merge annotations
annotation <- left_join(ann_locations, ann_genes, by = "cpgid")

# Merge with EWAS results
if (stratified == "no") {
    ewas <- left_join(ewas, annotation, by = "cpgid")
    ewas <- ewas[order(ewas$bacon.pval), ]
} else {
    ewas <- left_join(ewas, annotation, by = c("MarkerName" = "cpgid"))
    ewas <- ewas %>% dplyr::select(-Allele1, -Allele2)
    ewas$"P-value" <- as.numeric(ewas$"P-value")
    ewas <- ewas[order(ewas$"P-value"), ]
}

# Write annotated results
file_name <- paste0(out_dir, assoc, "_ewas_annotated_results", out_type)
fwrite(ewas, file = file_name)
