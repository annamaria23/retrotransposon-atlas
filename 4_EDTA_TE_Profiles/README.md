# General TE profiles of drosophilid species

The following scripts were used to generate TE profiles for 248 drosophilid species. 

To generate [RepeatModeler](https://github.com/Dfam-consortium/RepeatModeler) TE libraries:

``` bash

BuildDatabase -name ${BASE}/RM_db/genome_db ${BASE}/genome.fa
mkdir -p ${BASE}/RM_res && cd ${BASE}/RM_res
RepeatModeler -database ${BASE}/RM_db/genome_db -threads 4 -LTRStruct 


```
