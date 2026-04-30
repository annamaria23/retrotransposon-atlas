# Lifestyles of _Gypsy_-family transposons shape their regulatory mechanisms

Code repository for "Lifestyles of _Gypsy_-family transposons shape their regulatory mechanisms" by [Anna-Maria Papameletiou](https://github.com/annamaria23), Benjamin Czech Nicholson, [Susanne Bornelöv](https://github.com/susbo) and Gregory J Hannon. 

We present all the code used to annotate _D. melanogaster_-like _Gypsy_-family TEs across 248 drosophild genomes and downstream analysis.

<p align="center">
<img src="/Figures/Pipeline.png" width="500"/>
</p>

## Overview 

- Scripts to download and process the drosophilid genomes (from [Kim _et al_, 2021](https://elifesciences.org/articles/66405) and [Kim _et al_, 2024](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3002697)) are found [here](https://github.com/susbo/Drosophila_unistrand_clusters/tree/main/Genome_assemblies).
- _Drosophila melanogaster_ TE annotations : [Gypsy](1_Processing_Dmel_TEs/ALL_GYPSY_TES.fa) and all TEs
- Scripts to create hidden markov models (HMMs) for the POL, GAG, ENV and sORF2 open reading frames
- Scripts to generate TE profiles from drosophilid genomes using automated tools (RepeatMasker, EDTA).
- Scripts to annotate TEs
