# Load libraries
suppressPackageStartupMessages({
    library(R.utils)
    library(argparse)
    library(doFuture)
    library(progressr)
    library(progress)
    library(foreach)
    library(vars)
    library(dplyr)
    library(data.table)
    library(fst)
    library(tibble)
    library(tictoc)
    library(cli)
})

# Load functions
source("scripts/fxns/chunk_fxns.R")

############################################################################
# COMMAND LINE ARGUMENTS
############################################################################
parser <- argparse::ArgumentParser(description="Run EWAS comparing kidney tissue types")
parser$add_argument("--pheno", required=TRUE, help="Path to phenotype CSV or FST")
parser$add_argument("--methyl", required=TRUE, help="Path to methylation data CSV or FST")
parser$add_argument("--assoc", required=TRUE, help="Association variable (e.g. Type)")
parser$add_argument('--chunk-size', type="integer", default=1000, help="Number of CpGs per chunk")
parser$add_argument('--processing-type', '-pt', type="character", default="sequential", choices=c("sequential", "multisession", "multicore", "cluster"), help="Parallelization type")
parser$add_argument('--workers', type="integer", default=1, help="Number of parallel workers")
parser$add_argument('--out-dir', type="character", default="~/", help="Output directory")
parser$add_argument('--out-type', type="character", choices=c(".csv", ".csv.gz"), default=".csv", help="Output file type")
parser$add_argument('--out-prefix', type="character", default="all", help="Output prefix")

args <- parser$parse_args()

# Assign arguments
pheno_path <- args$pheno
mvals_path <- args$methyl
assoc_var <- args$assoc
chunk_size <- args$chunk_size
pt <- args$processing_type
n_workers <- args$workers
out_dir <- args$out_dir
out_type <- args$out_type
out_prefix <- args$out_prefix

threads_fst(nr_of_threads = n_workers)

############################################################################
# LOAD DATA
############################################################################
# Load phenotype
if (endsWith(pheno_path, ".fst")) {
    pheno <- read_fst(pheno_path)
} else {
    pheno <- fread(pheno_path)
}

# Ensure Type is treated as a factor
pheno$Type <- as.factor(pheno$Type)

# Load methylation data
if (endsWith(mvals_path, ".fst")) {
    mvals <- read_fst(mvals_path) %>% column_to_rownames(var=colnames(.)[1])
} else {
    mvals <- fread(mvals_path) %>% column_to_rownames(var=colnames(.)[1])
}

############################################################################
# CHECK INPUTS
############################################################################
if (!(assoc_var %in% colnames(pheno))) {
    stop(paste("Association variable", assoc_var, "not found in phenotype data."))
}

############################################################################
# CHUNK METHYLATION DATA
############################################################################
cpgs <- cpg.chunks(chunk_size, colnames(mvals))
mvals <- chunk.df(mvals, cpgs)

############################################################################
# EWAS FUNCTION 
############################################################################
registerDoFuture()
if (pt == "sequential") {
    plan(strategy = pt)
} else {
    plan(strategy = pt, workers = n_workers)
    options(future.globals.maxSize = +Inf)
}

ewas <- function(mvals, pheno) {
    p <- progressor(along = 1:length(mvals))
    results <- foreach(ii = 1:length(mvals), .combine = "rbind", .packages = c("vars")) %dopar% {
        m.chunk <- mvals[[ii]]
        p()
        foreach(i = colnames(m.chunk), .combine = "rbind", .packages = c("vars")) %dopar% {
            .GlobalEnv$m.chunk <- m.chunk
            # Model: CpG M-value ~ association variable (e.g., Type)
            formula_string <- paste0("m.chunk$", i, " ~ ", assoc_var)
            fit <- glm(formula = formula_string, data = pheno) %>%
                broom::tidy(conf.int = FALSE) %>%
                dplyr::filter(term == assoc_var) %>%
                dplyr::mutate(cpgid = i)
            fit
        }
    }
    return(results)
}

############################################################################
# RUN EWAS
############################################################################
handlers(global = TRUE)
results <- ewas(mvals, pheno)

print("EWAS complete.")
print(head(results))

############################################################################
# SAVE RESULTS
############################################################################
filename <- paste0(out_dir, out_prefix, "_", assoc_var, "_ewas_results", out_type)
if (endsWith(filename, ".gz")) {
    gz <- gzfile(filename, "wb")
    write.csv(results, gz, row.names = FALSE)
    close(gz)
} else {
    write.csv(results, file = filename, row.names = FALSE)
}

print(paste("EWAS results saved to:", filename))
