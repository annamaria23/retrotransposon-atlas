# Lifestyles of _Gypsy_-family transposons shape their regulatory mechanisms

Code repository for "Lifestyles of _Gypsy_-family transposons shape their regulatory mechanisms" by [Anna-Maria Papameletiou](https://github.com/annamaria23), Benjamin Czech Nicholson, [Susanne Bornelöv](https://github.com/susbo) and Gregory J Hannon. 

We present all the code used to annotate _D. melanogaster_-like _Gypsy_-family TEs across 248 drosophild genomes and downstream analysis.

## Overview 

### Data pre-procssing

<p align="center">
<img src="/Figures/figure1.png" width="600"/>
</p>

- Scripts to download and process the drosophilid genomes (from [Kim _et al_, 2021](https://elifesciences.org/articles/66405) and [Kim _et al_, 2024](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002697)) can be found [here](https://github.com/susbo/Drosophila_unistrand_clusters/tree/main/Genome_assemblies).
- Scripts to select genomes used in this study and generate statistics (N50, BUSCO completeness) can be found [here](./3_Genome_Selection_and_Preparation)
- Scripts to generate
 

<p align="center">
<img src="/Figures/Pipeline.png" width="500"/>
</p>

- **A-D**: Scripts to generate hidden markov models (HMMs) for the POL, GAG, ENV and sORF2 open reading frames can 



- _Drosophila melanogaster_ TE annotations : [Gypsy](1_Processing_Dmel_TEs/ALL_GYPSY_TES.fa) and all TEs
- Scripts to create hidden markov models (HMMs) for the POL, GAG, ENV and sORF2 open reading frames
- Scripts to generate TE profiles from drosophilid genomes using automated tools (RepeatMasker, EDTA).
- Scripts to annotate TEs
