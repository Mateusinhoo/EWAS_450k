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
if (!requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg38", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg38", ask = FALSE)
}

# Load 450k annotation from Bioconductor
suppressPackageStartupMessages({
  library(IlluminaHumanMethylation450kanno.ilmn12.hg38)
})

data("Locations")   # CpG positions
data("Other")       # UCSC and gene context
data("Manifest")    # Misc info

# Convert Locations to a data frame with cpgid column
ann_locations <- as.data.frame(Locations)
ann_locations$cpgid <- rownames(ann_locations)

# Optionally add UCSC gene info (e.g., gene names)
ann_genes <- as.data.frame(Other$UCSC_RefGene_Name)
colnames(ann_genes) <- "gene"
ann_genes$cpgid <- rownames(ann_genes)

# Merge into a final annotation table
annotation <- dplyr::left_join(ann_locations, ann_genes, by = "cpgid")

if(stratified=="no"){
    ewas <- left_join(ewas, annotation, by = "cpgid")
    ewas <- ewas[order(ewas$bacon.pval),]
} else{
        ewas <- left_join(ewas, annotation, by = c("MarkerName"="cpgid"))
        ewas <- ewas %>% dplyr::select(-Allele1, -Allele2)
        ewas$"P-value" <- as.numeric(ewas$"P-value")
        ewas <- ewas[order(ewas$"P-value"),]
}
file_name <- paste0(out_dir, assoc, "_ewas_annotated_results", out_type)
fwrite(ewas, file = file_name)
