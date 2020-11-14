#!/bin/R
library(MouseGastrulationData)
library(velociraptor)
library(scuttle)
library(uwot)
library(scater)

####
# Mixed gastrulation sample
####

md = AtlasSampleMetadata
samp = md[md$stage == "mixed_gastrulation", "sample"][1]

sce = EmbryoAtlasData(samp=samp, get.spliced = TRUE)
sce = sce[,!(sce$doublet|sce$stripped)]
sce = logNormCounts(sce)
reducedDim(sce, "new_umap") = uwot::umap(reducedDim(sce, "pca.corrected"), min_dist = 0.7)

velo = scvelo(
    sce,
    assay.X = "counts",
    assay.spliced = "spliced_counts",
    assay.unspliced = "unspliced_counts",
    use.dimred = "pca.corrected"
)

embedded <- embedVelocity(reducedDim(sce, "new_umap"), velo)
grid.df <- gridVectors(reducedDim(sce, "new_umap"), embedded)
sce$velocity_pseudotime <- velo$velocity_pseudotime

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

plotReducedDim(sce, "new_umap", colour_by="celltype") +
    geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, 
        xend=end.1, yend=end.2), arrow=arrow(length=unit(0.05, "inches"))) +
    scale_colour_manual(values = MouseGastrulationData::EmbryoCelltypeColours) +
    theme(legend.position = "none")

####
# Gut data
####
