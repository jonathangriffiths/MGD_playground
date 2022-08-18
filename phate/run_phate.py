#!/bin/python
# singularity shell ../images/phate.sif

import phate
import scanpy as sc
import pandas as pd
dat = sc.read_h5ad("../data/all/h5ad.h5ad")
model_default = phate.PHATE(n_pca=None)
phate_default = model_default.fit_transform(dat.obsm["pca.corrected"])
phate_pandas = pd.DataFrame(phate_default)
phate_pandas.index = dat.obs.index
phate_pandas.to_csv("phate_default_all.tsv", sep = "\t", header=False)
