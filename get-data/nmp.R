#!/bin/R
library(MouseGastrulationData)
library(ExperimentHub)
setExperimentHubOption("CACHE", "../cache")
library(LoomExperiment)
library(Matrix)
library(scuttle)

md = AtlasSampleMetadata
samp = md[md$stage == "E8.5", "sample"]

sce = EmbryoAtlasData(samp=samp, get.spliced = TRUE)
sce = sce[, sce$celltype %in% c("NMP", "Somitic mesoderm", 
    "Paraxial mesoderm", "Spinal cord", "Caudal Mesoderm", 
    "Caudal epiblast")]
sce = sce[,!(sce$doublet|sce$stripped)]
sce = logNormCounts(sce)

out = "../data/nmp"
if(!dir.exists(out)) dir.create(out, recursive=TRUE)
saveRDS(sce, file = file.path(out, "sce.rds"))

reducedDims(sce) = NULL #some bug in the reducedDims at the moment
scle <- as(sce, "SingleCellLoomExperiment")
LoomExperiment::export(scle, 
    file.path(out, "loom.loom"))

