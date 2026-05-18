#!/bin/bash
#SBATCH --cpus-per-task 1
#SBATCH --mem 2G
#SBATCH --time 0-01:00:00
#SBATCH --job-name orfs
#SBATCH --output log/%j-orfs.log

FILENAME=$(basename ${1})
OUTPUT=${2}

conda activate hmmer

#ENV
nhmmer --tblout ${OUTPUT}/${FILENAME}.env.hmm.tsv hmm/env.hmm ${1}
cat ${OUTPUT}/${FILENAME}.env.hmm.tsv  | grep -v "#" |  awk '$13 < 1e-10 &&  ((($8 - $7) < 0 ? -($8 - $7) : ($8 - $7)) > 100)' |  awk '{start=($7<$8 ? $7 : $8); end=($7<$8 ? $8 : $7); print $1, start-1, end, $3, $13, $12}' OFS="\t" > ${OUTPUT}/${FILENAME}.env.hmm.bed
bedtools getfasta -fi ${1} -bed ${OUTPUT}/${FILENAME}.env.hmm.bed -s >  ${OUTPUT}/${FILENAME}.env.hmm.fa
awk 'BEGIN{OFS="\t"} {print $1, "hmm_env", "CDS", $2+1, $3, ".", $5, "0", "gene_id \"" $4 "\";"}' ${OUTPUT}/${FILENAME}.env.hmm.bed > ${OUTPUT}/${FILENAME}.env.hmm.gtf

#GAG
nhmmer --tblout ${OUTPUT}/${FILENAME}.gag.hmm.tsv hmm/gag.hmm ${1}
cat ${OUTPUT}/${FILENAME}.gag.hmm.tsv | grep -v "#" | awk '$13 < 1e-20 && (($8 - $7 < 0 ? -($8 - $7) : $8 - $7) > 1000)' | awk '($1 !~ /^#/ && NF >= 13) {chr=$1; start=($7<$8 ? $7 : $8); end=($7<$8 ? $8 : $7); name=$3; strand=$12; print chr, start-1, end, name, strand}' OFS="\t" > ${OUTPUT}/${FILENAME}.gag.hmm.bed
bedtools getfasta -fi ${1} -bed ${OUTPUT}/${FILENAME}.gag.hmm.bed -s >  ${OUTPUT}/${FILENAME}.gag.hmm.fa
awk 'BEGIN{OFS="\t"} {print $1, "hmm_gag", "CDS", $2+1, $3, ".", $5, "0", "gene_id \"" $4 "\";"}' ${OUTPUT}/${FILENAME}.gag.hmm.bed > ${OUTPUT}/${FILENAME}.gag.hmm.gtf

#POL
nhmmer --tblout ${OUTPUT}/${FILENAME}.pol.hmm.tsv hmm/pol.hmm ${1}
cat ${OUTPUT}/${FILENAME}.pol.hmm.tsv | grep -v "#" | awk '$13 < 1e-20 && (($8 - $7 < 0 ? -($8 - $7) : $8 - $7) > 1000)' | awk '($1 !~ /^#/ && NF >= 13) {chr=$1; start=($7<$8 ? $7 : $8); end=($7<$8 ? $8 : $7); name=$3; strand=$12; print chr, start-1, end, name, strand}' OFS="\t" > ${OUTPUT}/${FILENAME}.pol.hmm.bed
bedtools getfasta -fi ${1} -bed ${OUTPUT}/${FILENAME}.pol.hmm.bed -s >  ${OUTPUT}/${FILENAME}.pol.hmm.fa
awk 'BEGIN{OFS="\t"} {print $1, "hmm_pol", "CDS", $2+1, $3, ".", $5, "0", "gene_id \"" $4 "\";"}' ${OUTPUT}/${FILENAME}.pol.hmm.bed > ${OUTPUT}/${FILENAME}.pol.hmm.gtf

# sORF2
blastx -query ${1} -subject env/brennecke_sorf2.fa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen" -max_target_seqs 1 -max_hsps 1 -evalue 1e-5 -num_threads 1 | awk '$3 > 50' | awk '$4 > 50' | awk '{ start = ($7 < $8 ? $7 : $8); end = ($7 > $8 ? $7 : $8); strand = ($7 < $8 ? "+" : "-"); print $1 "\tblastx_sorf2\tCDS\t" start "\t" end "\t.\t" strand "\t0\tgene_id \"" $2 "\";" }' > ${OUTPUT}/${FILENAME}.sorf2.gtf

# Full 
seqkit fx2tab -l -n ${1} | awk '{print $1 "\tfasta\tgene\t1\t" $2 "\t.\t+\t.\tgene_id \"" $1 "\";"}' > ${OUTPUT}/${FILENAME}.full.gtf

# LTRs
blastn -query ${1} -subject ${1} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen" \
| awk '$1 == $2' \
| awk '{
  qseqid=$1; qstart=$7; qend=$8; sstart=$9; send=$10;

  # normalize coordinates
  if (qstart > qend) { tmp=qstart; qstart=qend; qend=tmp }
  if (sstart > send) { tmp=sstart; sstart=send; send=tmp }

  key1=qstart":"qend":"sstart":"send;
  key2=sstart":"send":"qstart":"qend;

  if ( (qseqid SUBSEP key2) in seen ) {
    # retrieve previous hit
    split(seen[qseqid,key2], prev, " ");
    pqstart=prev[7]; pqend=prev[8];
    if (pqstart > pqend) { tmp=pqstart; pqstart=pqend; pqend=tmp }

    # ensure non-overlapping
    if ( (qend < pqstart) || (pqend < qstart) ) {

      # compute distance between starts (farther apart pair wanted)
      distance = (qstart < pqstart ? pqstart - qstart : qstart - pqstart);

      if (distance > bestdist[qseqid]) {
        bestdist[qseqid] = distance;
        best1[qseqid] = seen[qseqid,key2];
        best2[qseqid] = $0;
      }
    }
  } else {
    seen[qseqid,key1]=$0;
  }
}
END {
  for (seq in best1) {
    print best1[seq];
    print best2[seq];
  }
}' \
| awk '{print $1"\tLTR\tCDS\t" $7 "\t" $8 "\t.\t+\t0\tgene_id \"" $2 "\";"}' > ${OUTPUT}/${FILENAME}.LTRs.gtf
