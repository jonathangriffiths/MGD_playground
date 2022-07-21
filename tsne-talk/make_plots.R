#!/bin/R
library(MouseGastrulationData)
library(scran)
library(ggplot2)
library(irlba)
library(cowplot)
theme_set(theme_cowplot())
library(Rtsne)

sce = EmbryoAtlasData(samples=21)
sce = sce[,!sce$doublet & !sce$stripped]
sce = scuttle::logNormCounts(sce)

dec <- modelGeneVar(sce)
# plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
# curve(metadata(dec)$trend(x), col="blue", add=TRUE)
top.hvgs <- getTopHVGs(dec, prop=0.1)

pca = prcomp_irlba(t(logcounts(sce)[top.hvgs,]), n = 201)

## plot PCA

p = ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2, fill = sce$celltype)) +
    geom_point(pch = 21) +
    scale_fill_manual(values = MouseGastrulationData::EmbryoCelltypeColours) +
    theme(legend.position = "none")
ggsave(p, file = "pc1_vs_pc2.png", width = 4, height = 4, bg = "white")

p = ggplot(as.data.frame(pca$x), aes(x = PC3, y = PC4, fill = sce$celltype)) +
    geom_point(pch = 21) +
    scale_fill_manual(values = MouseGastrulationData::EmbryoCelltypeColours) +
    theme(legend.position = "none")
ggsave(p, file = "pc3_vs_pc4.png", width = 4, height = 4, bg = "white")

p = ggplot(as.data.frame(pca$x), aes(x = PC200, y = PC201, fill = sce$celltype)) +
    geom_point(pch = 21) +
    scale_fill_manual(values = MouseGastrulationData::EmbryoCelltypeColours) +
    theme(legend.position = "none")
ggsave(p, file = "pc200_vs_pc201.png", width = 4, height = 4, bg = "white")

## Plot distances
dm = as.matrix(dist(pca$x))
dim(dm)

cell = 42
dvec = dm[cell, seq_len(ncol(dm))[-cell]]
p = ggplot(mapping=aes(x=seq_along(dvec), y = dvec[order(dvec)])) +
    geom_bar(stat = "identity", fill = "black") +
    labs(x = "cells", y = "Distance to cell_42") +
    theme(axis.text.x = element_blank())
ggsave(p, file = "distance_to_cell42.png", width = 7, height = 4, bg = "white")

sdfac=30
simvec = exp(-(dvec^2)/sdfac)
p = ggplot(mapping=aes(x=seq_along(dvec), y = simvec[order(simvec, decreasing=TRUE)])) +
    geom_bar(stat = "identity", fill = "black") +
    labs(x = "cells", y = "exp(-d^2)") +
    theme(axis.text.x = element_blank())
ggsave(p, file = "-expd2_to_cell42.png", width = 7, height = 4, bg = "white")

xval = seq(from=0.01, to = 40, length.out = 500)
p=ggplot(mapping=aes(x=xval, y=exp(-xval^2/sdfac))) +
    geom_line() +
    labs(y = "'similarity' (exp(-d^2))", x = "distance")
ggsave(p, file = "dist_to_sim.png", width = 4, height = 4, bg = "white")

iter = 10^(0:4)
tsnes = lapply(iter, function(x) Rtsne(pca$x[, 1:50], pca=FALSE, max_iter=x))
plots = lapply(tsnes, function(x){
    ggplot(as.data.frame(x$Y), aes(x = V1, y =V2, fill = sce$celltype)) +
        geom_point(pch = 21) +
        scale_fill_manual(values = MouseGastrulationData::EmbryoCelltypeColours) +
        labs(x = 'tsne1', y = 'tsne2') +
        theme(legend.position="none") +
        ggtitle(paste0(x$max_iter, " iterations")) +
        theme(axis.ticks = element_blank(), axis.text=element_blank())

})
grid = plot_grid(plotlist = plots)
ggsave(grid, file = "tsne_grid_niter.png", bg = "white", width =12, height = 8)

perp = c(5, 10, 20, 50, 100, 200)
tsnes=lapply(perp, function(x) Rtsne(pca$x[, 1:50], pca=FALSE, perplexity=x))
plots = lapply(tsnes, function(x){
    ggplot(as.data.frame(x$Y), aes(x = V1, y =V2, fill = sce$celltype)) +
        geom_point(pch=21) +
        scale_fill_manual(values = MouseGastrulationData::EmbryoCelltypeColours) +
        labs(x = 'tsne1', y = 'tsne2') +
        theme(legend.position = "none") +
        ggtitle(paste0("perplexity=", x$perplexity)) +
        theme(axis.ticks = element_blank(), axis.text = element_blank())
})
grid = plot_grid(plotlist = plots)
ggsave(grid, file = "tsne_grid_perplexity.png", bg = "white", width = 12, height = 8)

p = ggplot(mapping=aes(x = sce$celltype, fill = sce$celltype)) +
    geom_bar(col = "black") +
    scale_fill_manual(values = MouseGastrulationData::EmbryoCelltypeColours) +
    theme(legend.position = "none",
        axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        axis.title.x = element_blank())
ggsave(p, file = "celltype_bar.png", bg = "white", width = 12, height = 5)

