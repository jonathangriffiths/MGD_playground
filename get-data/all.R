#!/bin/R
library(MouseGastrulationData)
library(ExperimentHub)
setExperimentHubOption("CACHE", "../cache")
library(Matrix)
library(scuttle)
library(zellkonverter)

sce = EmbryoAtlasData(get.spliced = TRUE)
sce = sce[, !(sce$doublet | sce$stripped)]
sce = logNormCounts(sce)

out = "../data/all"
if(!dir.exists(out)) dir.create(out, recursive=TRUE)
saveRDS(sce, file = file.path(out, "sce.rds"))

writeH5AD(sce, file = file.path(out, "h5ad.h5ad"), X_name = "logcounts")
#forget loom, it is doomed to the past
