#!/bin/bash

#SBATCH --cpus-per-task 1
#SBATCH --mem 1G
#SBATCH --time 0-00:30:00

conda activate minimap2 

BASEDIR=${1}
SPECIES=${2}
ASSEMBLY=${3}

BASE="${BASEDIR}/${SPECIES}/${ASSEMBLY}"

k8 /Users/papame01/calN50.js ${BASE}/genome.fa > ${BASE}/genome.calN50.stat

conda deactivate