import pandas as pd

# Load M-values (samples are columns, CpGs are rows)
mvals = pd.read_csv("data/mvals.csv.gz", index_col=0)

# Load phenotype files
pheno_blood = pd.read_csv("data/pheno_blood.csv")
pheno_brain = pd.read_csv("data/pheno_brain.csv")

# Get sample IDs
samples_blood = pheno_blood["sampleID"].tolist()
samples_brain = pheno_brain["sampleID"].tolist()

# Subset by columns (samples)
mvals_blood = mvals[samples_blood]
mvals_brain = mvals[samples_brain]

# Save to gzipped CSVs
mvals_blood.to_csv("data/mvals_blood.csv.gz", compression="gzip")
mvals_brain.to_csv("data/mvals_brain.csv.gz", compression="gzip")
