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
counts(sce) = as(counts(sce), "dgCMatrix")
sce = logNormCounts(sce)

out = "../data/nmp"
if(!dir.exists(out)) dir.create(out, recursive=TRUE)
saveRDS(sce, file = file.path(out, "sce.rds"))

loom = SingleCellLoomExperiment(sce)
LoomExperiment::export(loom, 
    file.path(out, "loom.loom"))
