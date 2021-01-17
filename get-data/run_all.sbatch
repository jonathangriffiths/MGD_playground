#!/bin/bash
#SBATCH -J get_mgd_data
#SBATCH -c 1
#SBATCH --mem 80G
#SBATCH -p main
#SBATCH -o slurm-%x-%j.out

files=($(ls *.R))
for file in ${files[@]}
do
    if [ ! -d "../data/$(basename $file .R)" ]
    then
        singularity exec --no-home -B .. ../images/mgd.sif Rscript $file
    fi
done
