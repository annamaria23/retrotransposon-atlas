# Genome selection and preparation 

Drosophilid were selected according to the following criteria: NN > 1,000,000 and N50 < 3,000. In cases where the same species had multiple genomes that passed the threshold, the highest quality was selected.

The [calN50](https://github.com/lh3/calN50) package was used to obtain AUC, N50, N80 and NN metrics for each genome:

```bash


k8 /Users/papame01/calN50.js ${BASE}/genome.fa > ${BASE}/genome.calN50.stat


```
