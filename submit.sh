#!/bin/bash -l

#$ -P me508
#$ -N g16
#$ -e g16.err
#$ -l h_rt=72:00:00
#$ -l mem_per_core=0.5G
#$ -pe omp 16
#$ -l cpu_arch=ivybridge

module load gaussian

g16 < 150K-ZERQOE-CO2.com  > 150K-ZERQOE-CO2.log
