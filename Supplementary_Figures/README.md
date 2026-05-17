# Data and scripts for the supplementary figures

This directory contains the scripts for creating, where relevant,the data files and scripts for creating the supplementary figures in the manuscript, summarised in the table below.

| **Figure** | **Data** | **Script** |
| :---: | :----: | :---: |
| **S1** | The file containing stats on the genomes used in **Figure S1A** can be found at [Species_Genome_Stats.csv](../3_Genome_Selection_and_Preparation/Species_Genome_Stats.csv) and the geographies of the species for **Figures S1B-C** were taken from [Pianezza _et al,_ 2025](https://github.com/rpianezza/Drosophilids-TE-biogeography/blob/main/D.melanogaster/tables/supplementary_file1.txt) |  R scripts for visualisation are at [tree_by_geog.R](./tree_by_geog.R)  |
| **S2** | Files containing the percentage of the genome covered by each TE superfamily per species can be found at [RM_EDTA_Results.csv](../4_EDTA_TE_Profiles/RM_EDTA_Results.csv) |  R scripts for visualisation are at [scatterplots.R](./scatterplots.R)  |
| **S3** | Files containing the percentage of the genome covered by each TE superfamily per species can be found at [RM_EDTA_Results.csv](../4_EDTA_TE_Profiles/RM_EDTA_Results.csv) and genome sizes can be found at [Species_Genome_Stats.csv](../3_Genome_Selection_and_Preparation/Species_Genome_Stats.csv) | The script for generating the phylogenetic tree can be found at [Phylogenetic_Trees1.R](../Main_Figures/Phylogenetic_Trees1.R) |
| **S5** | [LTR_similarity.txt](./LTR_similarity.txt) contains the file outputted as a result of blasting the front and back LTRs of a full-lenght TE against each other | [LTR_similarity_plot.R](./LTR_similarity_plot.R) contains the script to generate **Figure S5** |
|**S6** | Raw sRNA-seq reads were processed as described in the methods section and the mapping percentages can be found at [piRNA_mapping_percentages.csv](./piRNA_mapping_percentages.csv) | NA |
|**S7** | Raw sRNA-seq reads were processed as described in the methods section and a representative _.bam_ file was selected for each species for the rest of the analysis | [piRNA_lengths.R](piRNA_lengths.R) contains the script to create length plots and [piRNA_map.R](piRNA_map.R) contains the script to create logo plots for the anti-sense mapping piRNAs |
|**S8** | NA | NA |
|**S9** | NA | NA |
| **S10-S14** | NA | [TE_insertions_and_locations.R](TE_insertions_and_locations.R) contains the script to generate the plots relating to insertion counts of TEs in each species in **Figures S10-S14** |
|**S15** | NA | NA |
|**S16** | NA | NA |
|**S17-S18** | FIMO outputs can be found at (Main_Figures/all_TFs_LTRs_fimo.tsv)[../Main_Figures/all_TFs_LTRs_fimo.tsv] | Analysis scripts can be found at (Main_Figures/GLMM_with_TFs.R)[../Main_Figures/GLMM_with_TFs.R]  |
| **S19** | [ALL_sw_scores.tsv](ALL_sw_scores.tsv) contains a file with the SW scores of identified TFs in species compared to the ones in _D. melanogaster_| [TF_tree.R](TF_tree.R) contains the script to generate the phylogenetic tree in **Figure S19** |







