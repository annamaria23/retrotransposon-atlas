#!/bin/bash

conda activate mmseqs2

for FILE in /mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/250930_final/*/genome.fa; do
	
	SPECIES=$(echo "$FILE" | cut -d "/" -f 9)
	DIR=$(dirname "$FILE")

	cat ${DIR}/hmm_output.tsv | grep -v "#" | awk '$13 < 1e-20 && (($8 - $7 < 0 ? -($8 - $7) : $8 - $7) > 2000)' | awk '($1 !~ /^#/ && NF >= 13) {chr=$1; start=($7<$8 ? $7 : $8); end=($7<$8 ? $8 : $7); name=$3; strand=$12; print chr, start-1, end, name, strand}' OFS="\t" > ${DIR}/hmm_output.filt.bed

	# Get POL sequences
	bedtools getfasta -fi ${DIR}/genome.fa -bed ${DIR}/hmm_output.filt.bed -s >  ${DIR}/hmm_output.filt.bed.fa
	
	sbatch -p epyc -J ${SPECIES}_mmseqs --output logs/${SPECIES}_mmseqs.log 251013_mmseqs.sh ${DIR}

done
conda deactivate
