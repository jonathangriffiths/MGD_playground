#!/bin/python3

# Thanks to Valentine Svensson, from whom I "liberated" some code
# https://github.com/vals/Blog/blob/master/200917-trip-through-scvi/A%20trip%20through%20scVI.ipynb

####
# Module import
####

import pandas as pd
import anndata
import plotnine as p
import scvi
import scanpy as sc
import numpy as np

####
# Load and prepare data
####
fn = "../data/mixed_gast/loom.loom"
adata = anndata.read_loom(
    filename=fn
)
adata.var_names = adata.var.rownames
adata.obs_names = adata.obs.cell

####
# Fit scVI model
####

scvi.data.setup_anndata(
    adata,
    batch_key='sample'
)
model = scvi.model.SCVI(
    adata,
    gene_likelihood='nb',
    use_cuda=True
)
model.train(n_epochs=400)
model.save("mixed_gastrulation/")

####
# Visualise layers
####
latent = model.get_latent_representation()
adata.obsm["scvi_latent"] = latent

sc.pp.neighbors(adata, use_rep="scvi_latent")
sc.tl.umap(adata, min_dist=0.7)
sc.pl.umap(
    adata,
    color="celltype",
    frameon=False,
    save="mgast_latent.pdf"
)


####
# Check out denoised expression
####
sc.pl.umap(
    adata,
    color="T",
    layer="logcounts",
    gene_symbols="SYMBOL",
    frameon=False,
    save="mgast_T_logcounts.pdf"
)

sc.pl.umap(
    adata,
    color="T",
    gene_symbols="SYMBOL",
    frameon=False,
    save="mgast_T_counts.pdf"
)

adata.layers["scvi_normalised"] = model.get_normalized_expression(
    library_size=10e4
)
adata.layers["scvi_normalised_log"] = np.log2(adata.layers["scvi_normalised"] + 1)

sc.pl.umap(
    adata,
    layer = "scvi_normalised",
    color="T",
    gene_symbols="SYMBOL",
    frameon=False,
    save="mgast_T_denoised.pdf"
)

sc.pl.umap(
    adata,
    layer = "scvi_normalised_log",
    color="T",
    gene_symbols="SYMBOL",
    frameon=False,
    save="mgast_T_denoised_log2.pdf"
)



####
# Smoothed count histograms
####

g = 'T'
ens = adata.var[adata.var.SYMBOL == "T"].ENSEMBL[0]
adata.obs['T_counts']   = adata[:,ens].X.toarray()[:, 0]
adata.obs['T_logcounts']   = adata[:,ens].layers.get("logcounts").toarray()[:, 0]
adata.obs['T_smoothed']   = adata[:,ens].layers.get("scvi_normalised").toarray()[:, 0]
adata.obs['T_smoothed_lc']   = np.log2(adata.obs['T_smoothed'] + 1)

count_hist = adata.obs[f'{g}_counts'].value_counts().reset_index().rename(columns={'index': 'counts'})
p.options.figure_size = 6, 2
plot_ = (
    p.ggplot(p.aes(x='counts', y=f'{g}_counts'), count_hist.query('0 < counts < 25'))
    + p.geom_bar(stat='identity')
    + p.scale_x_log10()
    + p.theme_minimal()
    + p.labs(x=f'{g} UMI counts', y='Number cells')
)
plot_.save('mgast_T_counts.pdf', verbose=False)

count_hist = adata.obs[f'{g}_logcounts'].value_counts().reset_index().rename(columns={'index': 'logcounts'})
p.options.figure_size = 6, 2
plot_ = (
    p.ggplot(p.aes(x='logcounts'), count_hist.query('0 < logcounts < 25'))
    + p.geom_histogram(bins = 128, color = "k", fill = "w")
    + p.scale_x_log10()
    + p.theme_minimal()
    + p.labs(x=f'{g} UMI logcounts', y='Number cells')
)
plot_.save('mgast_T_logcounts.pdf', verbose=False)

count_hist = adata.obs[f'{g}_smoothed_lc'].value_counts().reset_index().rename(columns={'index': 'smoothed_lc'})
p.options.figure_size = 6, 2
plot_ = (
    p.ggplot(p.aes(x='smoothed_lc'), count_hist.query('0 < smoothed_lc < 25'))
    + p.geom_histogram(bins = 128, color = "k", fill = "w")
    + p.scale_x_log10()
    + p.theme_minimal()
    + p.labs(x=f'{g} UMI smoothed logcount', y='Number cells')
)
plot_.save('mgast_T_smoothed.pdf', verbose=False)


####
# Valentine's through-the-layers plots
###
# import torch
# from scvi.models.utils import one_hot
# px_conditional = []
# for tensors in full_posterior.sequential():
#     sample_batch, local_l_mean, local_l_var, batch_index, label = tensors
#     x_ = torch.log(1 + sample_batch)
#     qz_m, qz_v, z = full_posterior.model.z_encoder(x_)
#     ql_m, ql_v, library = full_posterior.model.l_encoder(x_)
#     px = full_posterior.model.decoder.px_decoder(qz_m, batch_index)
#     px_conditional.append(px.detach().cpu().numpy())
# px_conditional = np.vstack(px_conditional)

# tsne = TSNE(verbose=True, n_jobs=-1)
# YY = tsne.fit(px_conditional)
# for i, yy in enumerate(YY.T):
#     adata.obs[f'tsne_conditional_{i}'] = yy
