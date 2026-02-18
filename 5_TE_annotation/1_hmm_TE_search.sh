#!/bin/bash

conda activate hmmer

for FILE in /mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/250930_final/*/genome.fa; do
	
	SPECIES=$(echo "$FILE" | cut -d "/" -f 9)
	DIR=$(dirname "$FILE")
	
	sbatch -p epyc --cpus-per-task=1 --mem=24G --time=0-08:00:00 -J ${SPECIES}_hmmer_POL --output logs/${SPECIES}_hmmer_POL.log --wrap="nhmmer --tblout ${DIR}/hmm_output.tsv hmm/pol.hmm ${FILE}"

done
conda deactivate
