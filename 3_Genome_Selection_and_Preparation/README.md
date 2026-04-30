# Genome selection and preparation 

Drosophilid were selected according to the following criteria: NN > 1,000,000 and N50 < 3,000. In cases where the same species had multiple genomes that passed the threshold, the highest quality was selected.

The [calN50](https://github.com/lh3/calN50) tool was used to obtain AUC, N50, N80 and NN metrics for each genome:

```bash


k8 /Users/papame01/calN50.js ${BASE}/genome.fa > ${BASE}/genome.calN50.stat


```

[BUSCO](https://github.com/metashot/busco) completeness of the genomes using the _drosophilidae_odb12_ database was assessed by:


```bash

busco -i ${FILE} -l drosophilidae_odb12 -o busco_out --out_path ${DIR} -m genome --download_path /mnt/scratchc/ghlab/annamaria/important_files/busco_downloads --cpu 8 --offline


```

The resulting statistics for genomes used in this study can be found in [Species_Genome_Stats.csv](./Species_Genome_Stats.csv).

