# TE annotation pipeline

The following pipeline applies the _POL_ hidden Markov model as a starting point and outputs the consensus sequence for each _D. mel_-like _Gypsy_ TE in a specific species.


## Annotation of _Gypsy_-family TEs based on POL open reading frame sequences

Searching for all instances of POL in a genome:

```bash

nhmmer --tblout ${DIR}/hmm_output.tsv hmm/pol.hmm ${FILE}


```
Parsing the output file to select top POL hits:


```bash

cat ${DIR}/hmm_output.tsv | grep -v "#" | awk '$13 < 1e-20 && (($8 - $7 < 0 ? -($8 - $7) : $8 - $7) > 2000)' | awk '($1 !~ /^#/ && NF >= 13) {chr=$1; start=($7<$8 ? $7 : $8); end=($7<$8 ? $8 : $7); name=$3; strand=$12; print chr, start-1, end, name, strand}' OFS="\t" > ${DIR}/hmm_output.filt.bed # Format genomic locations

bedtools getfasta -fi ${DIR}/genome.fa -bed ${DIR}/hmm_output.filt.bed -s >  ${DIR}/hmm_output.filt.bed.fa # Extract genomic sequence


```

Lastly, the reciprocal best hit to the _D. melanogaster_ [reference](../1_Processing_Dmel_TEs/ALL_GYPSY_TES.fa) was obtained by running [run_mmseqs2.sh](./run_mmseqs2.sh). 


## Refinement of annotated _POL_ sequences into full consensus sequences

POL sequences were refined into full-length consensus sequences for each TE using a combination of tools, including [CIAlign](https://github.com/KatyBrown/CIAlign). The script for this is found in [3_make_consensus.sh](./3_make_consensus.sh). 


## Results

Results, including TE consensus sequence fasta files and downstream analyses can be found in [Resources](../Resources).
