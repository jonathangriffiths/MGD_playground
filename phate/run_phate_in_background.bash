#!/bin/bash
bind=$(dirname $(pwd))
image="../images/phate.sif"
nohup singularity exec -B $bind $image python run_phate.py &
