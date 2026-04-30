# RepeatModeller and EDTA runs on drosophilid genomes

## General TE profiles of drosophilid species

The following scripts were used to generate TE profiles for 248 drosophilid species. 

To generate [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) TE libraries:

``` bash

BuildDatabase -name ${BASE}/RM_db/genome_db ${BASE}/genome.fa
mkdir -p ${BASE}/RM_res && cd ${BASE}/RM_res
RepeatModeler -database ${BASE}/RM_db/genome_db -threads 4 -LTRStruct 


```

These libraries were then used to build [EDTA](https://github.com/oushujun/EDTA) libraries for the genomes:

``` bash

mkdir -p ${BASE}/EDTA_res && cd ${BASE}/EDTA_res

EDTA.pl --genome ${BASE}/genome.fa --step anno --anno 1 --threads 4 


```

EDTA was re-run on genomes for which contig naming prevented the porgramme from completing:


``` bash

mkdir -p ${BASE}/EDTA_res && cd ${BASE}/EDTA_res

if [[ ! -f "genome.fa.mod.EDTA.anno/genome.fa.mod.tbl" ]]; then

	seqkit replace -p '.*' -r 'seq{nr}' -o genome.renamed.fa genome.fa
	EDTA.pl --genome genome.renamed.fa --rmlib ../RM_res/*/consensi.fa.classified --species others --step all -- anno 1 --threads 4 --force 1 

fi


```

## RepeatModeller on _D. melanogaster_ TEs in each species 

``` bash

RepeatMasker -lib /mnt/scratchc/ghlab/annamaria/TE_annot_paper/250819_new_reference/ALL_GYPSY_TES.250826.fa -a -pa 1 -dir RM_Dmel/ genome.fa

```

## RepeatModller on species' curated TE libraries (to obtain putative insertion locations)

``` bash

RepeatMasker -lib 12_consensus/results_all_TEs.fa -a -pa 1 -dir RM_TElib/ genome.fa

```



