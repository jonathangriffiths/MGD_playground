#!/bin/R
library(MouseGastrulationData)
library(ExperimentHub)
setExperimentHubOption("CACHE", "../cache")
library(LoomExperiment)
library(Matrix)
library(scuttle)

sce = TChimeraData(samples = c(1:2, 5:10))
sce = sce[, !sce$celltype.mapped %in% c("Doublet", "Stripped")]
sce = logNormCounts(sce)

out = "../data/t-chim-85"
if(!dir.exists(out)) dir.create(out, recursive=TRUE)
saveRDS(sce, file = file.path(out, "sce.rds"))

reducedDims(sce) = NULL #some bug in the reducedDims at the moment
scle <- as(sce, "SingleCellLoomExperiment")
LoomExperiment::export(scle, 
    file.path(out, "loom.loom"))
