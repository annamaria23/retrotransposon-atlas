library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggtree)
library(treeio)
library(phytools)
library(ape)

full_sp <- read_delim("~/hannonlab Dropbox/Anna-Maria Papameletiou/TE_annot/SPECIES_NAME_CONVERSION.csv")
tree1 = read.tree("~/hannonlab Dropbox/Anna-Maria Papameletiou/TE_annot/Trees/gypsy10.fa.pol.hmm.fa.algn_cleaned.fasta.treefile") 
tree1$tip.label <- sub("^([^_]+)_.*", "\\1", tree1$tip.label)
tree1$tip.label <- full_sp$SpeciesName[match(tree1$tip.label, full_sp$Species)]
dupes <- tree1$tip.label[duplicated(tree1$tip.label)]
tree1 <- drop.tip(tree1, unique(dupes))

tree2= read.tree("~/hannonlab Dropbox/Anna-Maria Papameletiou/220705_phylogeny/input/All taxa_ML.con.tre") # Tree from Kim et al
dupes <- tree2$tip.label[duplicated(tree2$tip.label)]
tree2 <- drop.tip(tree2, unique(dupes))
common_tips <- intersect(tree1$tip.label, tree2$tip.label)

tree1 <- drop.tip(tree1, setdiff(tree1$tip.label, common_tips))
tree2 <- drop.tip(tree2, setdiff(tree2$tip.label, common_tips))

list <-  ENV %>% filter(TE == "gypsy10") %>% filter(!is.na(End_hmm_env)) %>%  filter(SpeciesName %in% tree1$tip.label) %>% pull(SpeciesName)

env_colours <- ifelse(tree1$tip.label %in% list, "magenta", "grey40")
cophy_obj <- cophylo(tree1, tree2, use.edge.length = F, rotate = T, tangle=tree1)

plot(cophy_obj, link.type = "curved",
     link.lty="solid",link.lwd=8,part = 0.46,
     link.col=make.transparent(env_colours,0.6),fzise=0.1, mar=c(0.1,0.1,2.1,0.1))
