import pandas as pd

# Load methylation matrix
mvals = pd.read_csv("data/mvals.csv.gz", index_col=0)

# Load phenotype files
pheno_blood = pd.read_csv("data/pheno_blood.csv")
pheno_brain = pd.read_csv("data/pheno_brain.csv")

# Get sample IDs
samples_blood = pheno_blood['sampleID'].tolist()
samples_brain = pheno_brain['sampleID'].tolist()

# Filter columns (samples) in mvals
mvals_blood = mvals[samples_blood]
mvals_brain = mvals[samples_brain]

# Save output
mvals_blood.to_csv("data/mvals_blood.csv.gz", compression='gzip')
mvals_brain.to_csv("data/mvals_brain.csv.gz", compression='gzip')

print(f"Saved blood: {mvals_blood.shape}")
print(f"Saved brain: {mvals_brain.shape}")
