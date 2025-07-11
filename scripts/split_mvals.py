import pandas as pd

# Load full M-values matrix (samples in rows, CpGs in columns)
mvals = pd.read_csv("data/mvals.csv.gz", index_col=0)

# Transpose to get CpGs as rows, samples as columns
mvals = mvals.T

# Load phenotype files
pheno_blood = pd.read_csv("data/pheno_blood.csv")
pheno_brain = pd.read_csv("data/pheno_brain.csv")

# Extract sample IDs
samples_blood = pheno_blood["sampleID"].tolist()
samples_brain = pheno_brain["sampleID"].tolist()

# Subset by sample ID (columns)
mvals_blood = mvals[samples_blood]
mvals_brain = mvals[samples_brain]

# Save the split methylation files
mvals_blood.to_csv("data/mvals_blood.csv.gz", compression="gzip")
mvals_brain.to_csv("data/mvals_brain.csv.gz", compression="gzip")
