#!/bin/R
# allowing homedir bind is important for accessing cache for Ehub/basilisk downloads
# singularity shell -B ~ ../images/mgd.sif
library(MouseGastrulationData)
library(zellkonverter)

meta = AtlasSampleMetadata
indices = sapply(split(meta, meta$stage), function(x) length(unique(x$pool_index)))
#5 E7.5 pools, the most
samples = meta$sample[meta$stage == "E7.5"]
sce = EmbryoAtlasData(samples = samples)
sce = sce[, !(sce$doublet | sce$stripped)]
sce = scuttle::logNormCounts(sce)
# quickly get hvgs
library(scran)
dec <- modelGeneVar(sce)
top.hvgs <- getTopHVGs(dec, prop=0.1)
rowData(sce)$is_hvg = rownames(sce)%in%top.hvgs

writeH5AD(sce, X_name="logcounts", file = "e75.h5ad")
