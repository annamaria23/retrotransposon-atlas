#!/bin/bash

#SBATCH --cpus-per-task 1
#SBATCH --mem 24G
#SBATCH --time 1-00:00:00

conda activate mmseqs2

mkdir -p ${1}/hmmer_db

# mmseqs2 run
mmseqs createdb ${1}/hmm_output.filt.bed.fa ${1}/hmmer_db/hmm_output_filt.db
mmseqs rbh --search-type 3 ${1}/hmmer_db/hmm_output_filt.db /mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/scripts/mmseqs_BL/BL_lib.fa.pol.nt.db ${1}/hmmer_db/hmm_output_pol_res.db ${1}/tmp
mmseqs convertalis ${1}/hmmer_db/hmm_output_filt.db /mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/scripts/mmseqs_BL/BL_lib.fa.pol.nt.db ${1}/hmmer_db/hmm_output_pol_res.db ${1}/hmm_pol_nt_output.tsv

#Adding output to file 
awk '!seen[$2]++' ${1}/hmm_pol_nt_output.tsv > ${1}/hmm_pol_nt_output.filt.tsv

# Turn into a bed file
cat ${1}/hmm_pol_nt_output.filt.tsv | awk '{
  split($1, a, "[:-]");
  chrom = a[1];
  start = a[2] - 1;
  end = a[3]; gsub(/\(\)/,"",end);
  name = $2;
  score = $12;
  strand = ".";
  print chrom"\t"start"\t"end"\t"name"\t"score"\t"strand 
}' > ${1}/hmm_pol_nt_output.filt.bed

# Get the sequence from the genome based on its location
bedtools getfasta -fi ${1}/genome.fa -bed ${1}/hmm_pol_nt_output.filt.bed -name -s > ${1}/hmm_pol_nt_output.filt.fa
