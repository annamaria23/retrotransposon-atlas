# Data and scripts for the main figures

This directory contains the scripts for creating and, where relevant,the data files for creating the figures in the manuscript, summarised in the table below.

| **Figure** | **Data** | **Script** |
| :---: | :----: | :---: |
| **1** | The data used in **Figure 1** is the output of the EDTA scripts and can be found at [RM_EDTA_Results.csv](../4_EDTA_TE_Profiles/RM_EDTA_Results.csv) | The script for generating phylogenetic trees can be found in  [Phylogenetic_Trees1.R](Phylogenetic_Trees1.R). The script for generating boxplots can be found at [boxplots.R](./boxplots.R) |
| **2** | The data used in **Figure 2** is the output of [...] | The script for generating phylogenetic trees can be found in  [TE_groups_plot.R](TE_groups_plot.R) |
| **3** | [...] | [...] |
| **4** | The data used to generate **Figure 4** is the [FIMO output for TE LTRs](all_TFs_LTRs_fimo.tsv). _ENV_ presense/absence was obtained from [ENV_locs.bed](ENV_locs.bed). To obtain a list of TFs expressed in the _D. melanogaster_ ovary, and specifically in somatic cells, germline cells or both, we obtained expression data from [FlyBase](https://wiki.flybase.org/wiki/FlyBase:Downloads_Overview#Single_Cell_RNA-Seq_Gene_Expression_.28scRNA-Seq_gene_expression_fb_.2A.tsv.gz.29).  | The script to run the GLMM model used in **Figure 4** and to visualise outputs can be found in [GLMM_with_TFs.R](GLMM_with_TFs.R) |




