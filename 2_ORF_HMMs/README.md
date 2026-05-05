## Obtaining HMM profiles from putative _Gypsy_ TEs

Input required:
- Predicted TE libraries for each genome assembly; we used the 101 drosophilid genomes (d101g): [B.Y. Kim et al., 2020](https://www.biorxiv.org/content/10.1101/2020.12.14.422775v1)
- GAG, POL and ENV peptide databases

#### 1. Obtaining TEs and processing them 

```bash

# Copying all TEs into one file
cat /mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/250516_automation_attempt/*/d101g/250626_all_consensus_seqs.fa | seqkit seq -m 3000 > all_d101g_TEs_v2.fa

# Assinging TE classes with TEsorter
conda activate tesorter
TEsorter all_d101g_TEs_v2.fa #The .lib file is the one of interest
conda deactivate

# Processing TEsorter output
# A random 1,000 TEs are selected
cat all_d101g_TEs_v2.fa.rexdb.cls.lib | seqkit grep -r -p "^.*#LTR/Gypsy$" | seqkit shuffle -s 42 | seqkit head -n 1000 | awk '/^>/{print $1; next}{print}' > TE_lib.fa

# blast against peptide databases

sbatch -p epyc --wrap "blastx -query TE_lib.fa -db ../../../../drosoTEs/blast_automation_attempt/scripts/env/pol.pep.new -max_target_seqs 1 -outfmt \"6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle\" > pol.blast.out"
sbatch -p epyc --wrap "blastx -query TE_lib.fa -db ../../../../drosoTEs/blast_automation_attempt/scripts/env/gag.pep.new -max_target_seqs 1 -outfmt \"6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle\" > gag.blast.out"
sbatch -p epyc --wrap "blastx -query TE_lib.fa -db ../../../../drosoTEs/blast_automation_attempt/scripts/env/env.pep -max_target_seqs 1 -outfmt \"6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle\" > env.blast.out"

# Extract sequences
bedtools getfasta -fi TE_lib.fa -bed <(awk '$3 > 40 && $4 > 300 {if ($5 < $6) print $1"\t"($5-1)"\t"$6"\t.\t.\t+"; else print $1"\t"($6-1)"\t"$5"\t.\t.\t-"}' pol.blast.out) -fo pol.fa
bedtools getfasta -fi TE_lib.fa -bed <(awk '$3 > 40 && $4 > 300 {if ($5 < $6) print $1"\t"($5-1)"\t"$6"\t.\t.\t+"; else print $1"\t"($6-1)"\t"$5"\t.\t.\t-"}' gag.blast.out) -fo gag.fa
bedtools getfasta -fi TE_lib.fa -bed <(awk '$3 > 40 && $4 > 300 {if ($5 < $6) print $1"\t"($5-1)"\t"$6"\t.\t.\t+"; else print $1"\t"($6-1)"\t"$5"\t.\t.\t-"}' env.blast.out) -fo env.fa

# MSA and HMM generation
conda activate hmmer
mafft pol.fa > pol.fa.aln && hmmbuild pol.hmm pol.fa.aln
mafft gag.fa > gag.fa.aln && hmmbuild gag.hmm gag.fa.aln
mafft env.fa > env.fa.aln && hmmbuild env.hmm env.fa.aln
conda deactivate

```

