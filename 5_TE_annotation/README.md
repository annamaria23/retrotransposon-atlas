# TE annotation pipeline

Searching for all instances of POL in a genome:

```bash

nhmmer --tblout ${DIR}/hmm_output.tsv hmm/pol.hmm ${FILE}


```
Parsing the output file to select top POL hits:


```bash

cat ${DIR}/hmm_output.tsv | grep -v "#" | awk '$13 < 1e-20 && (($8 - $7 < 0 ? -($8 - $7) : $8 - $7) > 2000)' | awk '($1 !~ /^#/ && NF >= 13) {chr=$1; start=($7<$8 ? $7 : $8); end=($7<$8 ? $8 : $7); name=$3; strand=$12; print chr, start-1, end, name, strand}' OFS="\t" > ${DIR}/hmm_output.filt.bed # Format genomic locations

bedtools getfasta -fi ${DIR}/genome.fa -bed ${DIR}/hmm_output.filt.bed -s >  ${DIR}/hmm_output.filt.bed.fa # Extract genomic sequence


```
