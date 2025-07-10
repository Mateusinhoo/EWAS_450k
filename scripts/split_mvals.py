import pandas as pd

# Load full methylation matrix (CpGs as rows, samples as columns)
mvals = pd.read_csv("data/mvals.csv.gz", index_col=0)

# Transpose so samples become rows, CpGs are columns
mvals = mvals.transpose()

# Load phenotypes
pheno_blood = pd.read_csv("data/pheno_blood.csv")
pheno_brain = pd.read_csv("data/pheno_brain.csv")

# Sample IDs from pheno files
blood_samples = pheno_blood["sampleID"].tolist()
brain_samples = pheno_brain["sampleID"].tolist()

# Subset rows (samples)
mvals_blood = mvals.loc[blood_samples]
mvals_brain = mvals.loc[brain_samples]

# Transpose back so CpGs are rows again and save
mvals_blood.transpose().to_csv("data/mvals_blood.csv.gz", compression="gzip")
mvals_brain.transpose().to_csv("data/mvals_brain.csv.gz", compression="gzip")

print("Saved blood methylation:", mvals_blood.shape, "→ data/mvals_blood.csv.gz")
print("Saved brain methylation:", mvals_brain.shape, "→ data/mvals_brain.csv.gz")
