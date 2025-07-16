# Make a BED file from EWAS results

suppressPackageStartupMessages({
    library(R.utils)
    library(argparse)
    library(data.table)
    library(tibble)
    library(dplyr)
})

# Arguments
parser <- argparse::ArgumentParser(description="Make BED file from EWAS results")
parser$add_argument("--results", required=TRUE, help="Path to EWAS annotated results file")
parser$add_argument("--out-dir", default="~/", help="Output directory")
parser$add_argument("--assoc", required=TRUE, help="Association variable used in EWAS")
args <- parser$parse_args()

# Load arguments
results <- args$results
out_dir <- args$out_dir
assoc <- args$assoc

# Load EWAS results
res <- fread(results)

# Check required columns
required_cols <- c("CpG_chrm", "CpG_beg", "CpG_end")
missing_cols <- setdiff(required_cols, colnames(res))
if (length(missing_cols) > 0) {
    stop(paste("ERROR: Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# Determine p-value column
if ("bacon.pval" %in% colnames(res)) {
    res <- res %>%
        select(CpG_chrm, CpG_beg, CpG_end, bacon.pval) %>%
        arrange(CpG_chrm, CpG_beg) %>%
        rename("#chrom" = CpG_chrm,
               start = CpG_beg,
               end = CpG_end,
               pvals = bacon.pval)
} else if ("P-value" %in% colnames(res)) {
    res <- res %>%
        select(CpG_chrm, CpG_beg, CpG_end, `P-value`) %>%
        arrange(CpG_chrm, CpG_beg) %>%
        rename("#chrom" = CpG_chrm,
               start = CpG_beg,
               end = CpG_end,
               pvals = `P-value`)
} else {
    stop("ERROR: No p-value column found (bacon.pval or P-value).")
}

# Save BED file
file_name <- paste0(out_dir, assoc, "_ewas_results.bed")
write.table(res, file = file_name, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

cat("BED file created at:", file_name, "\n")
