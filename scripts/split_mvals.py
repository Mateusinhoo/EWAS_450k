import pandas as pd

# Read brain phenotype and filter methylation data
pheno_brain = pd.read_csv("data/pheno_brain.csv")
mvals = pd.read_csv("data/mvals.csv.gz", index_col=0)
samples_brain = pheno_brain["sampleID"].tolist()
mvals_brain = mvals[samples_brain]
mvals_brain.to_csv("data/mvals_brain.csv.gz", compression="gzip")

# Read blood phenotype and filter methylation data
pheno_blood = pd.read_csv("data/pheno_blood.csv")
samples_blood = pheno_blood["sampleID"].tolist()
mvals_blood = mvals[samples_blood]
mvals_blood.to_csv("data/mvals_blood.csv.gz", compression="gzip")

print("Split complete. Saved to data/mvals_brain.csv.gz and data/mvals_blood.csv.gz")
