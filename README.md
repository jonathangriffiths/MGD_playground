# MGD_playground

This is a place for some simple code that demonstrates models etc. run on the MGD single-cell -omics data.

Why use _MouseGastrulationData_?
It has well annotated:

* Multi-sample, multi-timepoint scRNA-seq data from the gastrulation atlas
* Perturbation scRNA-seq data from the chimaeras
* (soon) scATAC-seq data from E8.25 embryos (with matched scRNA-seq in the atlas)
* (soon) Spatial transcriptomic data from E8.5 embryos.

all of which can be accessed using ExperimentHub and saved as loom files.

## Usage guide/structure

1. Make the singularity images from `images`

2. Prepare the data with `get-data` scripts

3. Run the scripts with the relevant image!
Scripts for a method are contained in their 
The naming should be self-explanatory (e.g. `bioc.sif` is for the R code). 

## Contributions

Please, if you have a method you would like to demonstrate, please do contribute it here - pull requests are very welcome.
Just make a folder that contains the scripts to run your code, and a definition file for a singularity image (or a docker image).