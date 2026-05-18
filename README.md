# Lifestyles of _Gypsy_-family transposons shape their regulatory mechanisms

Code repository for "Lifestyles of _Gypsy_-family transposons shape their regulatory mechanisms" by [Anna-Maria Papameletiou](https://github.com/annamaria23), Benjamin Czech Nicholson, [Susanne Bornelöv](https://github.com/susbo) and Gregory J Hannon. 

We present all the code used to annotate _D. melanogaster_-like _Gypsy_-family TEs across 248 drosophild genomes and associated downstream analysis.

## Overview 

### Data pre-procssing

<p align="center">
<img align="center" src="/Figures/figure1.png" width="550" hspace="20"/>
</p>

- Scripts to download and process the drosophilid genomes (from [Kim _et al_, 2021](https://elifesciences.org/articles/66405) and [Kim _et al_, 2024](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002697)) is the same as [van Lopik _et al_, 2023](https://www.nature.com/articles/s41467-023-42787-1) and can be found [here](https://github.com/susbo/Drosophila_unistrand_clusters/tree/main/Genome_assemblies).
- Scripts to select genomes used in this study and generate statistics (N50, BUSCO completeness) can be found [here](./3_Genome_Selection_and_Preparation). 
- Scripts to generate [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker) and [EDTA](https://github.com/oushujun/EDTA) profiles of the genomes can be found [here](./4_EDTA_TE_Profiles). 

 ### TE annotation pipeline

<p align="center">
<img src="/Figures/Pipeline.png" width="500"/>
</p>

- **A-D**: Scripts to generate hidden markov models (HMMs) for the POL, GAG, ENV and sORF2 open reading frames can be found [here](./2_ORF_HMMs).
- **E-H**: Scripts to identify _D. mel_-like _Gypsy_-family TEs in genomes and create high-quality consensus sequences can be found [here](./5_TE_annotation).

### Downstream analysis 

- Data and scripts used to generate the main figures in the paper can be found in [Main_Figures](Main_Figures).
- Data and scripts used to generate the spplementary figures in the paper can be found in [Supplementary_Figures](./Supplementary_Figures).


## Citing this work
...
