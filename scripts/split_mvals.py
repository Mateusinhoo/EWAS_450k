import pandas as pd

# Load the full M-values matrix (samples in rows, CpGs in columns)
mvals = pd.read_csv("data/mvals.csv.gz", index_col=0)

# Load phenotype files
pheno_blood = pd.read_csv("data/pheno_blood.csv")
pheno_brain = pd.read_csv("data/pheno_brain.csv")

# Extract sample IDs
samples_blood = pheno_blood["sampleID"].tolist()
samples_brain = pheno_brain["sampleID"].tolist()

# Subset rows (samples), then transpose back to CpGs as rows
mvals_blood = mvals.loc[samples_blood].T
mvals_brain = mvals.loc[samples_brain].T

# Save the split methylation files
mvals_blood.to_csv("data/mvals_blood.csv.gz", compression="gzip")
mvals_brain.to_csv("data/mvals_brain.csv.gz", compression="gzip")
