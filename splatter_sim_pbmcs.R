#### Loading packages
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scRNAseq)
  library(splatter)
  library(zellkonverter)
})

# Load h5ad file into a SingleCellExperiment object
sce <- readH5AD("../projects/Greta_clustering_errors/Data/pbmcs_filtered.h5ad",X_name = "counts")

# prepare counts for parameter estimation
counts_pbmc3k <- as.matrix(counts(sce))

# Estimate parameters from pbmc3k dataset
params <- splatEstimate(counts_pbmc3k, params = newSplatParams())
sim <- splatSimulate(params,group.prob = c(0.5, 0.5),method = "groups",verbose=FALSE)

# Save the SCE object to a .h5ad or .rds file
writeH5AD(sim, "../projects/Greta_clustering_errors/Data/sim_pbmcs_2groups.h5ad", X_name = "counts")
saveRDS(sim, file = "../projects/Greta_clustering_errors/Data/sim_pbmcs_2groups.rds")