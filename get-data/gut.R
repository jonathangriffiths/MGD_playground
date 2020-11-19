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
counts(sce) = as(counts(sce), "dgCMatrix")
sce = logNormCounts(sce)

out = "../data/gut"
if(!dir.exists(out)) dir.create(out, recursive=TRUE)
saveRDS(sce, file = file.path(out, "sce.rds"))

loom = SingleCellLoomExperiment(sce)
LoomExperiment::export(loom, 
    file.path(out, "loom.loom"))
