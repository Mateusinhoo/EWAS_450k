import pandas as pd

# Load full methylation matrix
mvals = pd.read_csv("data/mvals.csv.gz", index_col=0).transpose()

# Load full phenotype
pheno = pd.read_csv("data/pheno.csv")

# Split by cancer type
for tissue in ["blood", "brain"]:
    pheno_subset = pheno[pheno["cancer_type"] == tissue]
    sample_ids = pheno_subset["sampleID"].tolist()
    
    # Subset methylation matrix to matching samples
    mvals_subset = mvals.loc[sample_ids]
    
    # Transpose back so CpGs are rows and samples are columns
    mvals_subset = mvals_subset.transpose()
    
    # Save split methylation matrix and phenotype table
    mvals_subset.to_csv(f"data/mvals_{tissue}.csv.gz")
    pheno_subset.to_csv(f"data/pheno_{tissue}.csv", index=False)
