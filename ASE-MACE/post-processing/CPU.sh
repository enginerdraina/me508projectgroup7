#!/bin/bash -l
#$ -P me508 
#$ -N ASEMACE
#$ -e avcorr.err
#$ -l h_rt=12:00:00

module purge
module load miniconda
conda activate ASEMaceEnvironmentGPU

python oscillation_calcs.py

exit 

