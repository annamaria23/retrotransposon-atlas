library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggtree)
library(treeio)
library(paletteer)
library(ggtreeExtra)
library(forcats)
library(ggforce)
library(ggnewscale)
library(ape)

my_cols <- paletteer_d("ggthemes::Miller_Stone")
my_cols2 <- paletteer_d("ggsci::planetexpress_futurama")


tree <- read.iqtree("~/hannonlab Dropbox/Anna-Maria Papameletiou/trees_Kim2024/full_astral_fullAnnotation.tree")
tree@phylo$tip.label <- tree@phylo$tip.label %>%
  str_replace_all("_", ".") %>%
  str_to_sentence()
tips_to_drop <- setdiff(tree@phylo$tip.label, rm_sats$SpeciesName)
newtree <- drop.tip(tree@phylo, tips_to_drop)
p_tree <- ggtree(newtree, branch.length = "branch.length", layout = "circular") + theme_tree()

geog <- read_delim("~/hannonlab Dropbox/Anna-Maria Papameletiou/TE_annot/Supplement/rpianezza_supplementary_file1_geography.txt") %>% 
  setNames(c("assembly", "genus", "species", "SpeciesName", "realm1", "region")) %>% 
  mutate(region = ifelse(is.na(region), "Unknown", region))

region_colors <- setNames(
  c("#ea7580", "#f6a1a5", "#f8cd9c", "#1db6af", "#098bbe", "#172869", "black"),
  unique(geog$region)
)

clean_geog <- geog %>%
  mutate(region = ifelse(is.na(region), "Unknown", region))
region_list <- split(clean_geog$SpeciesName, clean_geog$region)
tree_grouped <- groupOTU(tree@phylo, region_list)
region_colors["0"] <- "grey85" 
if(!"Unknown" %in% names(region_colors)) {
  region_colors["Unknown"] <- "black"
}

ggtree(tree_grouped, aes(color = group), layout = "circular", size = 0.8) +
    scale_color_manual(values = region_colors, name = "Region") +
    geom_tippoint(aes(color = group), size = 1.5) +
  theme(legend.position = "right")
