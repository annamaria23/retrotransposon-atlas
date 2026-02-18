#!/bin/bash

#SBATCH --cpus-per-task 1
#SBATCH --mem 4G
#SBATCH --time 1-00:00:00
#SBATCH --job-name 12
#SBATCH --output log/%j.log

# Author: Susanne Bornelöv
# Last change: 2025-05-27 by Anna-Maria 

# Remap repeat models back to genome, extend region, cluster, and construct consensus sequences
# This is only applied for LTR transposons of high quality
# It works really (surprisingly!) well

# AMP: change to only add TEs that don't already exist 
# AMP: Change for new, POL-centric way to find TEs: blast POL agaisnt genome - extend 1,000 - consensus

conda activate phylo


mkdir -p ${1}/blast
mkdir -p ${1}/12_consensus
mkdir -p ${1}/12_consensus/out

# Blast database
makeblastdb -in ${1}/genome.fa -dbtype nucl -out ${1}/genome.fa

# Blast against genome and get top 20 best hits with 98% similarity
blastn -query ${1}/hmm_pol_nt_output.filt.fa -db ${1}/genome.fa -outfmt "6 qseqid sseqid pident length mismatch qstart qend sstart send sstrand qlen" -evalue 1e-20 -max_target_seqs 20 | awk '{OFS="\t"; if ($4 > $11/2) {print $0}}' | awk '$3 >= 90' > ${1}/12_consensus/out.blast.o

cat ${1}/12_consensus/out.blast.o | awk '{OFS="\t"; if ($1~/^\#/) {} else { if ($10~/plus/) {print $2, $8, $9, $1, $3, "+"} else {print $2, $9, $8, $1, $3, "-"}}}' | sort -k4,4 -k1,1 -k2,2n | uniq > ${1}/12_consensus/out.blast.bed

# Extend +/- 5,000 bp
bedtools slop -s -i ${1}/12_consensus/out.blast.bed -g ${1}/chrom.sizes -b 5000 > ${1}/12_consensus/out.blast.flank.bed


# Loop through identified TEs
TEs=`cat ${1}/hmm_pol_nt_output.filt.fa | grep ">" | sed 's/^.//' | awk -v OFS="" '{print $1}'`


for TE in $TEs
do
    outname=`echo $TE | cut -d':' -f1`
    outname="${outname}"
    mkdir -p ${1}/12_consensus/out/${outname}


	# Get sequences
    bedtools getfasta -fi ${1}/genome.fa -fo ${1}/12_consensus/out/${outname}/${outname}.blast.bed.fa -bed <(cat ${1}/12_consensus/out.blast.flank.bed | grep $TE) -s

    # Only if 3 or more entries (copies in genome) vecause othewise mafft doesn't work 
    count=$(grep -c '^>' "${1}/12_consensus/out/${outname}/${outname}.blast.bed.fa")

    if [ "$count" -gt 3 ]; then

	    # MSA
        mafft --reorder --auto --thread 1  ${1}/12_consensus/out/${outname}/${outname}.blast.bed.fa >  ${1}/12_consensus/out/${outname}/${outname}.maf.fa

	    # CIAlign
        CIAlign --infile ${1}/12_consensus/out/${outname}/${outname}.maf.fa --outfile_stem ${1}/12_consensus/out/${outname}/${outname}.maf.fa --crop_divergent --crop_divergent_min_prop_nongap 0.8 --crop_divergent_min_prop_ident 0.8 --remove_divergent --remove_divergent_minperc 0.3 --crop_ends --remove_insertions --insertion_max_size 3000 --remove_short --plot_input --plot_output --make_similarity_matrix_output --make_consensus

        # Remove flanking low-coverage regions from consensus≈
        blastn -query ${1}/12_consensus/out/${outname}/${outname}.maf.fa_consensus.fasta -db ${1}/genome.fa -outfmt 6 > ${1}/12_consensus/out/${outname}/out.blast.o

	    #Remove if just making a log file -- check if this is ok
	    perl /mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/scripts/Susanne/data/crop_zero_coverage.pl ${1}/12_consensus/out/${outname}/out.blast.o ${1}/12_consensus/out/$outname/$outname.maf.fa_consensus.fasta > ${1}/12_consensus/out/$outname/$outname.cons.log

    else
        seqkit head -n 1 ${1}/12_consensus/out/${outname}/${outname}.blast.bed.fa > ${1}/12_consensus/out/${outname}/${outname}.cons.fa 
    fi


    HEADER="${2}_${outname}"
    sed "s/^>.*/>${HEADER}/" ${1}/12_consensus/out/${outname}/${outname}.cons.fa > ${1}/12_consensus/out/${outname}/${outname}.cons.fixname.fa

    bash 251014_get_LTRs.sh ${1}/12_consensus/out/${outname}/${outname}.cons.fixname.fa ${1}/12_consensus/out/${outname}/${outname}.cons.LTRs.gtf
    awk '{print $1"\t" $4-1 "\t" $5 "\t" $1"_TE\t.\t+"}' ${1}/12_consensus/out/${outname}/${outname}.cons.LTRs.gtf > ${1}/12_consensus/out/${outname}/${outname}.cons.LTRs.bed

    if [ -s "${1}/12_consensus/out/${outname}/${outname}.cons.LTRs.bed" ]; then

        bedtools getfasta -fi ${1}/12_consensus/out/${outname}/${outname}.cons.fixname.fa -bed ${1}/12_consensus/out/${outname}/${outname}.cons.LTRs.bed -name -s |  sed "s/^>.*/>${HEADER}/" > ${1}/12_consensus/out/${outname}/${outname}.final.sequence.fa

    else

        scp ${1}/12_consensus/out/${outname}/${outname}.cons.fixname.fa ${1}/12_consensus/out/${outname}/${outname}.final.sequence.fa

    fi

done

conda deactivate