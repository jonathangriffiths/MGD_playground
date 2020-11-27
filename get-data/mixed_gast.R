#!/bin/R
library(MouseGastrulationData)
library(ExperimentHub)
setExperimentHubOption("CACHE", "../cache")
library(LoomExperiment)
library(Matrix)
library(scuttle)

md = AtlasSampleMetadata
samp = md[md$stage == "mixed_gastrulation", "sample"][1]

sce = EmbryoAtlasData(samp=samp, get.spliced = TRUE)
sce = sce[,!(sce$doublet|sce$stripped)]
for(ass in names(assays(sce))){
    assay(sce, ass) = as(assay(sce, ass), "dgCMatrix")
}
# counts(sce) = as(counts(sce), "dgCMatrix")
sce = logNormCounts(sce)

out = "../data/mixed_gast"
if(!dir.exists(out)) dir.create(out, recursive=TRUE)
saveRDS(sce, file = file.path(out, "sce.rds"))

sce2 = sce
reducedDims(sce2) = NULL
scle <- as(sce2, "SingleCellLoomExperiment")
LoomExperiment::export(scle, 
    file.path(out, "loom.loom"))
