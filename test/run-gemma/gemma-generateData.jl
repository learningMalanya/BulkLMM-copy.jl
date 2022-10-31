# Will be using BXD data
## Include the script to generate data
include("../../src/readData.jl");

## Read in BXD data:
pheno_file = "../../data/bxdData/BXDtraits.csv"
pheno = readBXDpheno(pheno_file);
geno_file = "../../data/bxdData/BXDgeno_prob.csv"
geno = readGenoProb_ExcludeComplements(geno_file);


## Read in BXD data to the BIMBAM format:
transform_bxd_geno_to_gemma(geno_file, "BXDgeno_prob_bimbam.txt");
transform_bxd_pheno_to_gemma(pheno_file, "BXDpheno_bimbam.txt", 1997);