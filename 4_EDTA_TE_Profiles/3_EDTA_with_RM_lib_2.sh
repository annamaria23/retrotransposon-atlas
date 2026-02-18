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

if [[ ! -f "genome.fa.mod.EDTA.anno/genome.fa.mod.tbl" ]]; then

	seqkit replace -p '.*' -r 'seq{nr}' -o genome.renamed.fa genome.fa
	EDTA.pl --genome genome.renamed.fa --rmlib ../RM_res/*/consensi.fa.classified --species others --step all -- anno 1 --threads 4 --force 1 

fi



conda deactivate
