import pandas as pd

# Read the full methylation matrix (samples might be in rows, CpGs in columns)
mvals = pd.read_csv("data/mvals.csv.gz", index_col=0)

# Transpose so samples are columns
mvals = mvals.T

# Split brain
pheno_brain = pd.read_csv("data/pheno_brain.csv")
samples_brain = pheno_brain["sampleID"].tolist()
mvals_brain = mvals[samples_brain].T
mvals_brain.to_csv("data/mvals_brain.csv.gz", compression="gzip")

# Split blood
pheno_blood = pd.read_csv("data/pheno_blood.csv")
samples_blood = pheno_blood["sampleID"].tolist()
mvals_blood = mvals[samples_blood].T
mvals_blood.to_csv("data/mvals_blood.csv.gz", compression="gzip")

print(" Done splitting methylation file into brain and blood.")
