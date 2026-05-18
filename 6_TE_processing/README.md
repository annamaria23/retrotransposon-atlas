# TE sequences processing: phylogenetic tree building and gtf file creation

## Obtaining gtf files with ORF locations

TE gtf files were obtained using the following script:

```bash

TEs=$(cat ../results/BL_POL.names.txt)

conda activate phylo

mkdir -p ../results/by_TE_FINAL

for te in ${TEs}; do
	echo ${te}	
	mkdir -p ../results/by_TE_FINAL/${te}

	sbatch -p epyc -J ${te}_phylo -o logs_teannot/${te}_phylo.out -e logs_teannot/${te}_phylo.err --cpus-per-task=1 --mem=2G --time=0-02:00:00 --wrap="
	source /Users/papame01/miniconda3/etc/profile.d/conda.sh;
	# General extraction of ORFs
	conda activate hmmer;
	seqkit grep -r -p "_${te}\$" ../results/by_TE_FINAL/annotated_TEs.fa > ../results/by_TE_FINAL/${te}/${te}.fa;
	bash getgtfs.sh /mnt/scratche/slow/ghlab/annamaria/drosoTEs/blast_automation_attempt/results/by_TE_FINAL/${te}/${te}.fa /mnt/scratche/slow/ghlab/annamaria/drosoTEs/blast_automation_attempt/results/by_TE_FINAL/${te};
	conda deactivate;

	#Tree building based on POL
	conda activate phylo;
	cd ../results/by_TE_FINAL/${te}/;
	mafft ${te}.fa.pol.hmm.fa > ${te}.fa.pol.hmm.fa.algn;
	CIAlign --infile  ${te}.fa.pol.hmm.fa.algn --outfile_stem ${te}.fa.pol.hmm.fa.algn --crop_divergent --crop_divergent_min_prop_nongap 0.8 --crop_divergent_min_prop_ident 0.8 --remove_divergent --remove_divergent_minperc 0.3 --crop_ends --remove_insertions --insertion_max_size 3000 --remove_short --plot_input --plot_output --make_similarity_matrix_output --make_consensus;
	iqtree -s ${te}.fa.pol.hmm.fa.algn_cleaned.fasta -m MFP -bb 1000 -alrt 1000 -abayes -B 1000 -T AUTO -bnni --ninit 100 --seed 42 -redo;
	conda deactivate;

"
done


```


## Data preparation for dN/dS calculations

TE ORFs were processed for dN/dS calculations for the _POL_, _ENV_, and _sORF2_ open reading frames using the following script:

```bash

conda activate phylo

FILENAME="$(basename "$1")"
FILEPATH="$(dirname "$1")"
ORF=${2}

# minsize of ORF, depending on ORF we are using
case "${ORF}" in
    "POL")
        MINSIZE=500
        ;;
    "ENV")
        MINSIZE=100
        ;;
    "sORF2")
        MINSIZE=50
        ;;
esac

echo ${ORF}
echo ${MINSIZE}

mkdir -p ${FILEPATH}/${FILENAME}.dnds

awk '/^>/ {sub(/:.*/, "", $0)} {print}' ${1} > ${FILEPATH}/${FILENAME}.dnds/${ORF}.rename.fa
getorf -sequence ${FILEPATH}/${FILENAME}.dnds/${ORF}.rename.fa -outseq ${FILEPATH}/${FILENAME}.dnds/${ORF}.nt.fa -find 3 -minsize ${MINSIZE}

awk '/^>/ {sub(/ .*/, "", $0)} {print}' ${FILEPATH}/${FILENAME}.dnds/${ORF}.nt.fa > ${FILEPATH}/${FILENAME}.dnds/${ORF}.nt.final.fa
seqkit rmdup -n ${FILEPATH}/${FILENAME}.dnds/${ORF}.nt.final.fa > ${FILEPATH}/${FILENAME}.dnds/${ORF}.nt.final.final.fa
transeq -frame 1 -sequence ${FILEPATH}/${FILENAME}.dnds/${ORF}.nt.final.final.fa -outseq ${FILEPATH}/${FILENAME}.dnds/aa.fa >/dev/null
mafft --thread 16 ${FILEPATH}/${FILENAME}.dnds/aa.fa > ${FILEPATH}/${FILENAME}.dnds/aa.fa.maf 2>/dev/null
pal2nal.pl ${FILEPATH}/${FILENAME}.dnds/aa.fa.maf ${FILEPATH}/${FILENAME}.dnds/${ORF}.nt.final.final.fa -output paml > ${FILEPATH}/${FILENAME}.dnds/aa.maf.phy
iqtree2 -s ${FILEPATH}/${FILENAME}.dnds/aa.fa.maf -m MFP -redo 

cd ${FILEPATH}/${FILENAME}.dnds
cp *.treefile ${FILENAME}tree.nwk

```
