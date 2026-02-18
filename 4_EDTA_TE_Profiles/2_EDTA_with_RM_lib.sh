#!/bin/bash

#SBATCH --cpus-per-task 4
#SBATCH --mem 80G
#SBATCH --time 10-00:00:00

conda activate EDTA 

unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
unset PERL_MB_OPT
unset PERL_MM_OPT

BASEDIR=${1}
SPECIES=${2}
ASSEMBLY=${3}

BASE="${BASEDIR}/${SPECIES}/${ASSEMBLY}"

mkdir -p ${BASE}/EDTA_res

cd ${BASE}/EDTA_res # I think this might be needed to specify the output path

EDTA.pl --genome ${BASE}/genome.fa --step anno --anno 1 --threads 4 

conda deactivate
