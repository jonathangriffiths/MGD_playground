#!/bin/R
library(MouseGastrulationData)
library(ExperimentHub)
setExperimentHubOption("CACHE", "../cache")
library(LoomExperiment)
library(Matrix)
library(scuttle)

sce = EmbryoAtlasData(get.spliced = TRUE)
sce = sce[, sce$celltype %in% c("Gut", "Visceral endoderm", 
    "ExE endoderm", "Def. endoderm", "Notochord")]
sce = sce[, !(sce$doublet | sce$stripped)]
sce = logNormCounts(sce)

out = "../data/gut"
if(!dir.exists(out)) dir.create(out, recursive=TRUE)
saveRDS(sce, file = file.path(out, "sce.rds"))

reducedDims(sce) = NULL #some bug in the reducedDims at the moment
scle <- as(sce, "SingleCellLoomExperiment")
LoomExperiment::export(scle, 
    file.path(out, "loom.loom"))

