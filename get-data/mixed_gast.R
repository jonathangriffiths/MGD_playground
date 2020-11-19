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
counts(sce) = as(counts(sce), "dgCMatrix")
sce = logNormCounts(sce)

out = "../data/mixed_gast"
if(!dir.exists(out)) dir.create(out, recursive=TRUE)
saveRDS(sce, file = file.path(out, "sce.rds"))

loom = SingleCellLoomExperiment(sce)
LoomExperiment::export(loom, 
    file.path(out, "loom.loom"))
