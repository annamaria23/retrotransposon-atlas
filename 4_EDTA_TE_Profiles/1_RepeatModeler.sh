#!/bin/bash

#SBATCH --cpus-per-task 4
#SBATCH --mem 80G
#SBATCH --time 6-00:00:00

conda activate repeat_modeler 

unset PERL5LIB
unset PERL_LOCAL_LIB_ROOT
unset PERL_MB_OPT
unset PERL_MM_OPT

BASEDIR=${1}
SPECIES=${2}
ASSEMBLY=${3}

BASE="${BASEDIR}/${SPECIES}/${ASSEMBLY}"

mkdir -p ${BASE}/RM_db
mkdir -p ${BASE}/RM_res


# Building genome database
BuildDatabase -name ${BASE}/RM_db/genome_db ${BASE}/genome.fa

cd ${BASE}/RM_res # I think this might be needed to specify the output path

RepeatModeler -database ${BASE}/RM_db/genome_db -threads 4 -LTRStruct 

conda deactivate