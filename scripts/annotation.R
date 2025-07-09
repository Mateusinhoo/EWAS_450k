# Script for running adding annotation data to ewas results
suppressPackageStartupMessages({
    library(argparse)
    library(dplyr)
    library(data.table)
})

# Define command line arguments
parser <- argparse::ArgumentParser(description="Script for adding annotation data to ewas results")
parser$add_argument('--input-file', '-i', required=TRUE, help="Path to ewas results data")
parser$add_argument('--out-dir', required=TRUE, help="Path to output directory")
parser$add_argument('--stratified', choices=c("yes", "no"), default="no", help="Results from a stratified analysis: yes or no")
parser$add_argument('--assoc', required=TRUE, type="character", help="Association variable EWAS was performed with.")
parser$add_argument('--out-type', type="character", choices=c(".csv", ".csv.gz"), default=".csv", help="Output file type")
parser$add_argument('--annotation', required=TRUE, help="Path to manifest annotation file")
parser$add_argument('--snp-key', required=TRUE, help="Path to SNP annotation file")
parser$add_argument('--out-prefix', type="character", default="all", help="Prefix for output files")

# Parse arguments
args <- parser$parse_args()
out_prefix <- args$out_prefix
results <- args$`input_file`
out_dir <- args$`out_dir`
stratified <- args$`stratified`
assoc <- args$`assoc`
out_type <- args$`out_type`
annotation_path <- args$`annotation`
snp_key_path <- args$`snp_key`

# Load files
ewas <- fread(results)
ann_locations <- fread(annotation_path)
ann_genes <- fread(snp_key_path)

# Assume CpG IDs are in "probeID", fallback if needed
if ("Probe_ID" %in% colnames(ann_genes)) {
  ann_genes <- ann_genes %>% rename(cpgid = Probe_ID)
} else {
  stop("No valid CpG ID column found in snp_annotation file.")
}

# Rename CpG ID column in ann_locations (manifest file)
if ("probeID" %in% colnames(ann_locations)) {
  ann_locations <- ann_locations %>% dplyr::rename(cpgid = probeID)
} else if ("IlmnID" %in% colnames(ann_locations)) {
  ann_locations <- ann_locations %>% dplyr::rename(cpgid = IlmnID)
} else {
  print(colnames(ann_locations))  # For debugging if it fails again
  stop("No valid CpG ID column found in annotation_manifest file.")
}

# Merge annotation data
annotation <- left_join(ann_locations, ann_genes, by = "cpgid")

# Merge annotation with EWAS results
if (stratified == "no") {
    ewas <- left_join(ewas, annotation, by = "cpgid")
    ewas <- ewas[order(ewas$bacon.pval), ]
} else {
    ewas <- left_join(ewas, annotation, by = c("MarkerName" = "cpgid"))
    ewas <- ewas %>% dplyr::select(-Allele1, -Allele2)
    ewas$`P-value` <- as.numeric(ewas$`P-value`)
    ewas <- ewas[order(ewas$`P-value`), ]
}

# Save output
filename <- paste0(out_dir, out_prefix, "_", assoc, "_ewas_annotated_results", out_type)
fwrite(ewas, file = filename)
