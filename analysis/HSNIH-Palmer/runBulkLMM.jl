using Helium;
using DataFrames;
using CSV;
using GeneNetworkAPI;

pheno_10441 = get_pheno("HSNIH-Palmer", "10441");
geno_filename = "HSNIH-Palmer.he";


geno_mat = readhe(geno_filename);
NaN_to_remove = zeros(size(geno_mat, 1));
for p = 1:size(geno_mat, 1)
    
    if any(!isfinite, geno_mat[p, :])
        
        NaN_to_remove[p] = 1;
        
    end
    
end

geno_mat = geno_mat[NaN_to_remove .!= 1, :];

match_id = Array{Int, 1}(undef, 6147);
for i in 1:6147
    
    indicator = 0;
    for j in 1:4099
        if geno_sampleID[i] == pheno_sampleID[j]
            indicator = indicator + 1;
        end
    end
    
    match_id[i] = indicator;
end
samples_in_common = geno_sampleID[match_id .==1];

match_id_pheno = Array{Int, 1}(undef, 4099);
for i in 1:length(match_id_pheno)
    
    indicator = 0;
    for j in 1:length(samples_in_common)
        if pheno_sampleID[i] == samples_in_common[j]
            indicator = indicator + 1;
        end
    end
    
    match_id_pheno[i] = indicator;
end

## subset geno:
geno_common = geno_mat[:, match_id .== 1];
## subset pheno:
pheno_common = pheno_10441[match_id_pheno .== 1, :];

n = size(geno_common, 2);
mean_freq = sum(geno_common; dims = 2)./(2*n);
less_common_markers = map(x -> (x > 0.90) || (x < 0.10), mean_freq)[:, 1];

geno_processed = permutedims(geno_common[less_common_markers .!= 1, :]);
pheno_processed = reshape(pheno_common[:, 4], :, 1);

kinship_Gc = CSV.read("output/kinship.cXX.txt", DataFrame, delim = '\t', header = false);
kinshipMat_Gc = Matrix(kinship_Gc); 

@time results_BulkLMM = scan(pheno_processed, geno_processed, kinshipMat_Gc; reml = false, method = "alt")