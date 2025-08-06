library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Input args from Snakemake
sample_sheet_path <- snakemake@input[["sample_sheet"]]
idat_dir <- snakemake@input[["idat_dir"]]
output_file <- snakemake@output[["mvals"]]

# Set working directory to IDATs
setwd(idat_dir)

# Load sample sheet
sample_sheet <- read.csv(sample_sheet_path)
targets <- sample_sheet[, c("Sentrix ID", "Sentrix Position")]
colnames(targets) <- c("Slide", "Array")
targets$Sample_Name <- sample_sheet$ID

# Read IDAT files
rgSet <- read.metharray.exp(targets = targets)

# Preprocess using Funnorm
mSet <- preprocessFunnorm(rgSet)

# Extract M-values
mvals <- getM(mSet)

# Save
write.csv(mvals, file = output_file)
