# Script for performing Epigenome-Wide association analysis that takes in command line arguments

# Import libraries
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

# Source functions from outside scripts
source("scripts/fxns/chunk_fxns.R")

################################################################################################################################
#                                   DEFINE AND PARSE COMMAND LINE ARGUMENTS                                                    # 
################################################################################################################################
# Define command line arguments
parser <- argparse::ArgumentParser(description="Script for running EWAS")
parser$add_argument("--pheno", required=TRUE, help="Path to phenotype data")
parser$add_argument("--methyl", required=TRUE, help="Path to methylation data")
parser$add_argument("--assoc", required=TRUE, type="character", help="Association variable")
parser$add_argument('--stratified', choices=c("yes", "no"), default="no", help="Stratified analysis: yes or no")
parser$add_argument("--chunk-size", type="integer", default=1000, help="Number of CpGs per chunk")
parser$add_argument('--processing-type', '-pt', type="character", default="sequential", choices=c("sequential", "multisession", "multicore", "cluster"), help="Parallelization type")
parser$add_argument('--workers', type="integer", default=1, help="Number of parallel workers")
parser$add_argument('--out-dir', type="character", default="~/", help="Output directory")
parser$add_argument('--out-type', type="character", choices=c(".csv", ".csv.gz"), default=".csv", help="Output file type")
parser$add_argument('--out-prefix', type="character", default="all", help="Output prefix")
parser$add_argument('--subset-condition', type="character", default=NULL, help="Subset condition like 'cancer_type == \"brain\"'")

# Parse arguments
args <- parser$parse_args()
pheno_path <- args$pheno
mvals_path <- args$methyl
assoc_var <- args$assoc
stratified <- args$stratified
chunk_size <- args$chunk_size
pt <- args$processing_type
n_workers <- args$workers
out_dir <- args$out_dir
out_type <- args$out_type
out_prefix <- args$out_prefix
subset_condition <- args$`subset-condition`

# Set number of threads for fst
threads_fst(nr_of_threads = n_workers)

####################################################################################################
#                                   READ IN DATA                                                   #
####################################################################################################
# Read phenotype
if (endsWith(pheno_path, ".fst")) {
  pheno <- read_fst(pheno_path) %>% column_to_rownames(var=colnames(.)[1])
} else {
  pheno <- fread(pheno_path) %>% column_to_rownames(var=colnames(.)[1])
}

# Apply subset condition if provided
if (!is.null(subset_condition)) {
  message("Subsetting phenotype data using: ", subset_condition)
  pheno <- subset(pheno, eval(parse(text = subset_condition)))
}

# Create sample_subgroup if missing
if (!"sample_subgroup" %in% colnames(pheno)) {
  pheno$sample_subgroup <- paste0(pheno$cancer_type, "_", pheno$ethnicity)
}
if (assoc_var == "sample_subgroup") {
  pheno$sample_subgroup <- factor(pheno$sample_subgroup)
}

# Read methylation
if (endsWith(mvals_path, ".fst")) {
  mvals <- read_fst(mvals_path) %>% column_to_rownames(var=colnames(.)[1])
} else {
  mvals <- fread(mvals_path) %>% column_to_rownames(var=colnames(.)[1])
}

####################################################################################################
#                                     CHECK DATA                                                   #
####################################################################################################
if (!(assoc_var %in% colnames(pheno))) {
  stop(paste("Association variable", assoc_var, "not found in phenotype data."))
} else {
  pheno <- pheno %>% relocate(all_of(assoc_var))
}

####################################################################################################
#                                      CHUNK DATA                                                  #
####################################################################################################
cpgs <- cpg.chunks(chunk_size, colnames(mvals))
mvals <- chunk.df(mvals, cpgs)

####################################################################################################
#                              LINEAR REGRESSION ANALYSIS                                          #
####################################################################################################
registerDoFuture()
if (pt == "sequential") {
  plan(strategy = pt)
  cat("Processing run sequentially.\n")
} else {
  plan(strategy = pt, workers = n_workers)
  options(future.globals.maxSize = +Inf)
  cat("Asynchronous parallel processing using", pt, "with", n_workers, "worker(s).\n")
}

cat("Starting EWAS: ", out_prefix, "\n")
tic()
handlers(global = TRUE)
if (interactive()) {
  handlers(list(handler_progress(format = "[:bar] :percent ELAPSED::elapsed, ETA::eta")))
} else {
  handlers("cli")
  options(cli.progress_handlers = "progressr")
}

ewas <- function(mvals, pheno) {
  tic()
  p <- progressor(along = 1:length(mvals))
  results <- foreach(ii = 1:length(mvals), .combine = "rbind", .packages = c("vars")) %dopar% {
    m.chunk <- mvals[[ii]]
    p()
    foreach(i = colnames(m.chunk), .combine = "rbind", .packages = c("vars")) %dopar% {
      .GlobalEnv$m.chunk <- m.chunk
      covariates <- c(assoc_var, "age", "smoking_status", "bmi")
      covariates <- covariates[covariates %in% colnames(pheno)]
      if (length(covariates) == 0) stop("None of the covariates (assoc_var, age, smoking_status, bmi) found in phenotype data.")
      model.base <- paste(covariates, collapse = " + ")
      string.formula <- paste0("m.chunk$", i, " ~ ", model.base)
      fit <- glm(formula = string.formula, data = pheno) %>%
        broom::tidy(conf.int = FALSE) %>%
        dplyr::filter(grepl(paste0("^", assoc_var), term)) %>%
        dplyr::mutate(cpgid = i)
      fit
    }
  }
  toc()
  return(results)
}

results <- ewas(mvals, pheno)

print("Number of EWAS results:")
print(nrow(results))
print("First few rows of results:")
print(head(results))

# Save results
filename <- paste0(out_dir, out_prefix, "_", assoc_var, "_ewas_results", out_type)
if (endsWith(filename, ".gz")) {
  gz <- gzfile(filename, "wb")
  write.csv(results, gz, row.names = FALSE)
  close(gz)
} else {
  write.csv(results, file = filename, row.names = FALSE)
}
