#!/bin/R
# singularity shell ../images/bioc.sif
library(MouseGastrulationData)
library(SingleCellExperiment)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

plotdir = "plots"
dir.create(plotdir, showWarnings=FALSE)

sce = readRDS("../data/all/sce.rds")
#only need the cell metadata plus the phate vis
cd = colData(sce)
rm(sce); gc()
phate = read.delim("phate_default_all.tsv", header=FALSE, row.names=1)

names(phate) = paste0("phate", 1:2)

cd = as.data.frame(cbind(cd, phate))

p = ggplot(cd[sample(nrow(cd)), ], aes(x = phate1, y = phate2, col = celltype)) +
    geom_point() +
    scale_colour_manual(values = EmbryoCelltypeColours) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "none")
ggsave(p, file = file.path(plotdir, "phate_default_allcells.png"), width = 8,
    height = 8, bg = "white")
#that kinda sucks
# maybe it is capturing global distances better (e.g. the big changes in blood/ExE endoderm,
# as well as Emb. vs. ExE) but it mangles the embryonic data terribly.
# seems like a greater focus on local structure provides something more usable.
