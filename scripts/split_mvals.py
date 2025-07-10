import pandas as pd

# Load phenotypes
pheno_blood = pd.read_csv("data/pheno_blood.csv")
pheno_brain = pd.read_csv("data/pheno_brain.csv")

# Load mvals normally (samples are in rows)
mvals = pd.read_csv("data/mvals.csv.gz")

# Filter to blood and brain samples
mvals_blood = mvals[mvals["sampleID"].isin(pheno_blood["sampleID"])]
mvals_brain = mvals[mvals["sampleID"].isin(pheno_brain["sampleID"])]

# Save each split
mvals_blood.to_csv("data/mvals_blood.csv.gz", index=False)
mvals_brain.to_csv("data/mvals_brain.csv.gz", index=False)
