library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggtree)
library(treeio)
library(paletteer)
library(ggtreeExtra)
library(forcats)
#library(ggforce)
library(ggnewscale)
library(ape)

mdg1_list <-  c("412", "blood", "mdg1", "Tabor", "Stalker", "Stalker4")
gypsy_list <- c("gypsy5", "ZAM", "Tirant", "accord", "accord2", "Idefix", "Quasimodo", "McClintock", "HMS-Beagle", "rover", "17.6", "297",
                "transpac", "chouto", "gypsy4", "gypsy10", "gypsy9", "Burdock", "gypsy7", "gypsy", "gtwin", "gypsy2", "gypsy6", "springer",
                "gypsy3", "HMS", "opus", "HMS-Beagle2", "gypsy1", "Transpac", "Shellder")
mdg3_list <- c("mdg3", "micropia", "invader", "Spoink", "invader3", "invader2", "invader1", "invader6", "MLE", "Bica")
osvaldo_list <- c("gypsy8", "gypsy12", "Circe", "Souslik")
chimpo_list <- c("Chimpo")

species <- read_delim("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/species_list.txt") %>% 
    mutate(fullname_1 = substr(`Full name`, 1, 1)) %>% 
    mutate(fullname_2 = sub("^[^ ]+ ", "", `Full name`)) %>% 
    mutate(SpeciesName = paste0(fullname_1, ".", fullname_2)) %>% 
    select(c("Species", "SpeciesName")) %>% 
    distinct()

more_species <- read_delim("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/298_species_list_NEW.txt", col_names = FALSE) %>% 
    select(c("X1", "X2")) %>% 
    setNames(c("Species", "SpeciesName")) %>% 
    separate(SpeciesName, into = c("Part1", "Part2"), sep = "_") %>% 
    mutate(
        Part1 = str_to_title(Part1),
        Part2 = tolower(Part2)
    ) %>% 
    mutate(SpeciesName = paste(Part1, Part2, sep = ".")) %>% 
    select(-c("Part1", "Part2"))

species_names <- rbind(species, more_species)

TEs_kept <- read_delim("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/results/251021_by_POL/annot_final.tsv", "\t", col_names = c("Name", "ORF", "Start", "End", ".", "pvalue", "0", "gene_id", "gene_name")) %>%
  mutate(Species = sub("_(?!.*_).*", "", Name, perl = TRUE)) %>% 
  mutate(TE = sub(".*_", "", Name)) %>% 
  left_join(species_names) %>%
  na.omit()%>%
  pull(SpeciesName)

TEs_mmseqs <- read_delim("ALL_species.txt", "\t", col_names = FALSE) %>%
  select(c("X2", "X3", "X13")) %>% 
  setNames(c("TE", "Sim", "Species")) %>% 
  mutate(Species = map_chr(Species, ~ str_split(.x, "/")[[1]][9])) %>% 
  mutate(TE = sub(":.*", "", TE)) %>% 
  filter(!TE %in% c("MLE", "gypsy3", "Stalker2", "Shellder")) %>%
  add_count(TE, name = "te_count") %>% 
  mutate(TE = fct_reorder(TE, te_count, .desc = TRUE)) %>% 
  group_by(TE, Species) %>%
  filter(Sim == max(Sim)) %>%
  slice(1) %>%
  ungroup()  %>%
  complete(Species, TE) %>% 
  left_join(species_names) 

TEs_mdg1 <- TEs_mmseqs %>% filter(TE %in% mdg1_list) %>% mutate(TE = fct_drop(TE))
TEs_gypsy <- TEs_mmseqs %>% filter(TE %in% gypsy_list) %>% mutate(TE = fct_drop(TE))
TEs_mdg3 <- TEs_mmseqs %>% filter(TE %in% mdg3_list) %>% mutate(TE = fct_drop(TE))
TEs_osvaldo <- TEs_mmseqs %>% filter(TE %in% osvaldo_list) %>% mutate(TE = fct_drop(TE))
TEs_chimpo <- TEs_mmseqs %>% filter(TE %in% chimpo_list) %>% mutate(TE = fct_drop(TE))

TE_counts_mdg1 <- table(TEs_mdg1 %>% group_by(TE, Species) %>% unique() %>% na.omit() %>% pull(SpeciesName)) %>% as.data.frame() %>% setNames(c("SpeciesName", "Count")) 
TE_counts_gypsy <- table(TEs_gypsy %>% group_by(TE, Species) %>% unique() %>% na.omit() %>% pull(SpeciesName)) %>% as.data.frame() %>% setNames(c("SpeciesName", "Count")) 
TE_counts_mdg3 <- table(TEs_mdg3 %>% group_by(TE, Species) %>% unique() %>% na.omit() %>% pull(SpeciesName)) %>% as.data.frame() %>% setNames(c("SpeciesName", "Count")) 
TE_counts_osvaldo <- table(TEs_osvaldo %>% group_by(TE, Species) %>% unique() %>% na.omit() %>% pull(SpeciesName)) %>% as.data.frame() %>% setNames(c("SpeciesName", "Count")) 
TE_counts_chimpo <- table(TEs_chimpo %>% group_by(TE, Species) %>% unique() %>% na.omit() %>% pull(SpeciesName)) %>% as.data.frame() %>% setNames(c("SpeciesName", "Count")) 


# Tree
#tree <- read.iqtree("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/Kim_tree.tree")
tree <- read.iqtree("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/full_astral_fullAnnotation.tree")
labels <- tree@phylo$tip.label
labels <- tolower(labels)
labels <- gsub("_", ".", labels)
substr(labels, 1, 1) <- toupper(substr(labels, 1, 1))
tree@phylo$tip.label <- labels

tips_to_drop <- c(setdiff(tree@phylo$tip.label, TEs_mmseqs$SpeciesName), "D.pachea")
newtree <- drop.tip(tree@phylo, tips_to_drop)
p_tree <- ggtree(newtree, branch.length = "branch.length", layout = "fan", open.angle=90) + theme_tree()
p_tree <- rotate_tree(p_tree, 90)

#Scale so same for all plots
scale_mdg1 <- scale_fill_gradient(low = "lightyellow", high = "#5e3300", na.value = "beige", limits = c(0.55, 1), name = "% Similarity")
scale_gypsy <- scale_fill_gradient(low = "lightyellow", high = "#005193", na.value = "beige", limits = c(0.55, 1), name = "% Similarity")
scale_mdg3 <- scale_fill_gradient(low = "lightyellow", high = "darkred", na.value = "beige", limits = c(0.55, 1), name = "% Similarity")
scale_osvaldo <- scale_fill_gradient(low = "lightyellow", high = "#125b00", na.value = "beige", limits = c(0.55, 1), name = "% Similarity")
scale_chimpo <- scale_fill_gradient(low = "lightyellow", high = "#43006c", na.value = "beige", limits = c(0.55, 1), name = "% Similarity")

p_tree +
  geom_fruit(
    data = TEs_mdg1,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.2, pwidth = 2,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1))+
    scale_color_manual(values = c("TRUE" = "#5e3300", "FALSE" = "transparent"), guide = "none") +
  scale_mdg1 + 
  new_scale_fill() + 
  geom_fruit(
    data = TE_counts_mdg1,
    geom = geom_col,
    orientation = "y",
    mapping = aes(y = SpeciesName, x = Count),
    fill = "#5e3300", alpha=0.5,
    offset = 0.3, pwidth = 1,
    axis.params = list(axis = "x", text.size = 3, text.angle = 0))
ggsave("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/260330_TE_mdg1_tree.pdf", device = cairo_pdf, height=6, width=6)

p_tree$data <- p_tree$data %>% mutate(dmel = label == "D.melanogaster")

p_tree +
  geom_tippoint(mapping = aes(subset = dmel), shape = 21, fill = "red", color = "black", size = 3, alpha = 1)+
  geom_fruit(
    data = TEs_gypsy,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.1, pwidth = 4.5,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1))+
  scale_gypsy + 
  new_scale_fill() + 
  scale_color_manual(values = c("TRUE" = "#005193", "FALSE" = "transparent"), guide = "none") +
  geom_fruit(
    data = TE_counts_gypsy,
    geom = geom_col,
    orientation = "y",
    mapping = aes(y = SpeciesName, x = Count),
    fill = "#005193", alpha=0.5,
    offset = 0.3, pwidth = 2,
    axis.params = list(axis = "x", text.size = 3, text.angle = 0))+
  geom_hilight(node=getMRCA(as.phylo(p_tree), c(c("D.ficusphila", "D.yakuba"))), fill="#638B66FF", alpha=0.7, extendto=1)+
  geom_hilight(node=getMRCA(as.phylo(p_tree), c(c("Z.bogoriensis", "Z.camerounensis"))), fill="#D7CE9FFF", alpha=0.7, extendto=1)
ggsave("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/260330_TE_gypsy_tree.pdf", device = cairo_pdf, height=14, width=14)

p_tree +
  geom_fruit(
    data = TEs_mdg3,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.2, pwidth = 2.3,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1))+
  scale_mdg3 + 
  scale_color_manual(values = c("TRUE" = "darkred", "FALSE" = "transparent"), guide = "none") +
  new_scale_fill() + 
  geom_fruit(
    data = TE_counts_mdg3,
    geom = geom_col,
    orientation = "y",
    mapping = aes(y = SpeciesName, x = Count),
    fill = "darkred", alpha=0.5,
    offset = 0.3, pwidth = 1,
    axis.params = list(axis = "x", text.size = 3, text.angle = 0))
ggsave("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/260330_TE_mdg3_tree.pdf", device = cairo_pdf, height=6.5, width=6.5)

p_tree +
  geom_fruit(
    data = TEs_osvaldo,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.3, pwidth = 1,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1))+
  scale_osvaldo + 
  scale_color_manual(values = c("TRUE" = "#125b00", "FALSE" = "transparent"), guide = "none") +
  new_scale_fill() + 
  geom_fruit(
    data = TE_counts_osvaldo,
    geom = geom_col,
    orientation = "y",
    mapping = aes(y = SpeciesName, x = Count),
    fill = "#125b00", alpha=0.5,
    offset = 0.4, pwidth = 0.8,
    axis.params = list(axis = "x", text.size = 3, text.angle = 0))
ggsave("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/260330_TE_osvaldo_tree.pdf", device = cairo_pdf, height=5, width=5)


p_tree +
  geom_fruit(
    data = TEs_chimpo,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.1, pwidth = 6,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1))+
  scale_chimpo + 
  scale_color_manual(values = c("TRUE" = "#43006c", "FALSE" = "transparent"), guide = "none") +
  new_scale_fill() + 
  geom_fruit(
    data = TE_counts_chimpo,
    geom = geom_col,
    orientation = "y",
    mapping = aes(y = SpeciesName, x = Count),
    fill = "#43006c", alpha=0.5,
    offset = 0.2, pwidth = 0.2,
    axis.params = list(axis = "x", text.size = 3, text.angle = 0))
ggsave("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/260330_TE_chimpo_tree.pdf", device = cairo_pdf, height=4, width=4)


subtree <- extract.clade(newtree, node = getMRCA(newtree, c("D.ficusphila", "D.yakuba")))
subtree$edge.length[is.na(subtree$edge.length)] <- 0
p_tree_sub <- ggtree(subtree, branch.length = "branch.length", layout = "rectangular") + theme_tree2() +  geom_tiplab(aes(label = ""))

p_tree_sub +
geom_tiplab(size = 4, align = TRUE) + 
    geom_fruit(
    data = TEs_gypsy,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 1.2, pwidth = 5,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_gypsy + 
  scale_color_manual(values = c("TRUE" = "#005193", "FALSE" = "transparent"), guide = "none") +
  new_scale_fill() + 
  new_scale_color()+
geom_fruit(
    data = TEs_mdg1,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.2, pwidth = 1.3,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_mdg1 + 
  scale_color_manual(values = c("TRUE" = "#5e3300", "FALSE" = "transparent"), guide = "none") +
  new_scale_fill() + 
  new_scale_color()+
  geom_fruit(
    data = TEs_mdg3,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.2, pwidth = 1.7,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_mdg3 + 
  scale_color_manual(values = c("TRUE" = "darkred", "FALSE" = "transparent"), guide = "none") +
 new_scale_fill() + 
 new_scale_color()+
  geom_fruit(
    data = TEs_osvaldo,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.2, pwidth = 0.7,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_osvaldo +
  scale_color_manual(values = c("TRUE" = "#125b00", "FALSE" = "transparent"), guide = "none") +
  new_scale_fill() + 
  new_scale_color()+
   geom_fruit(
    data = TEs_chimpo,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.2, pwidth = 7,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_chimpo+
  scale_color_manual(values = c("TRUE" = "#43006c", "FALSE" = "transparent"), guide = "none")
  ggsave("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/260330_TE_tree_mel_2.pdf", height=8, width=13, device = cairo_pdf)



subtree <- extract.clade(newtree, node = getMRCA(newtree, c("Z.bogoriensis", "Z.camerounensis")))
subtree$edge.length[is.na(subtree$edge.length)] <- 0
p_tree_sub <- ggtree(subtree, branch.length = "branch.length", layout = "rectangular") + theme_tree2() +  geom_tiplab(aes(label = ""))

p_tree_sub +
geom_tiplab(size = 4, align = TRUE) + 
    geom_fruit(
    data = TEs_gypsy,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 1.2, pwidth = 5,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_gypsy + 
  scale_color_manual(values = c("TRUE" = "#005193", "FALSE" = "transparent"), guide = "none") +
  new_scale_fill() + 
  new_scale_color()+
geom_fruit(
    data = TEs_mdg1,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.2, pwidth = 1.3,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_mdg1 + 
  scale_color_manual(values = c("TRUE" = "#5e3300", "FALSE" = "transparent"), guide = "none") +
  new_scale_fill() + 
  new_scale_color()+
  geom_fruit(
    data = TEs_mdg3,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.2, pwidth = 1.7,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_mdg3 + 
  scale_color_manual(values = c("TRUE" = "darkred", "FALSE" = "transparent"), guide = "none") +
 new_scale_fill() + 
 new_scale_color()+
  geom_fruit(
    data = TEs_osvaldo,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.2, pwidth = 0.7,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_osvaldo +
  scale_color_manual(values = c("TRUE" = "#125b00", "FALSE" = "transparent"), guide = "none") +
  new_scale_fill() + 
  new_scale_color()+
   geom_fruit(
    data = TEs_chimpo,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.2, pwidth = 7,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_chimpo+
  scale_color_manual(values = c("TRUE" = "#43006c", "FALSE" = "transparent"), guide = "none")
ggsave("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/260330_TE_tree_zap_2.pdf", height=8, width=13, device = cairo_pdf)


#tree <- read.iqtree("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/Kim_tree.tree")
tree <- read.iqtree("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/full_astral_fullAnnotation.tree")
labels <- tree@phylo$tip.label
labels <- tolower(labels)
labels <- gsub("_", ".", labels)
substr(labels, 1, 1) <- toupper(substr(labels, 1, 1))
tree@phylo$tip.label <- labels
tips_to_drop <- c(setdiff(tree@phylo$tip.label, TEs_mmseqs$SpeciesName), "D.pachea")
newtree <- drop.tip(tree@phylo, tips_to_drop)
p_tree <- ggtree(newtree, branch.length = "branch.length", layout = "rectangular") + theme_tree()

p_tree +
geom_tiplab(size = 4, align = TRUE) + 
    geom_fruit(
    data = TEs_gypsy,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 1.2, pwidth = 5,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_gypsy + 
  scale_color_manual(values = c("TRUE" = "#005193", "FALSE" = "transparent"), guide = "none") +
  new_scale_fill() + 
  new_scale_color()+
geom_fruit(
    data = TEs_mdg1,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.2, pwidth = 1.3,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_mdg1 + 
  scale_color_manual(values = c("TRUE" = "#5e3300", "FALSE" = "transparent"), guide = "none") +
  new_scale_fill() + 
  new_scale_color()+
  geom_fruit(
    data = TEs_mdg3,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.2, pwidth = 1.7,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_mdg3 + 
  scale_color_manual(values = c("TRUE" = "darkred", "FALSE" = "transparent"), guide = "none") +
 new_scale_fill() + 
 new_scale_color()+
  geom_fruit(
    data = TEs_osvaldo,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.2, pwidth = 0.7,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_osvaldo +
  scale_color_manual(values = c("TRUE" = "#125b00", "FALSE" = "transparent"), guide = "none") +
  new_scale_fill() + 
  new_scale_color()+
   geom_fruit(
    data = TEs_chimpo,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TE, fill = as.numeric(Sim), color=!is.na(Sim)),
    offset = 0.2, pwidth = 7,
    axis.params = list(axis = "x", text.size = 4 , text.angle = 90, hjust=1))+
  scale_chimpo+
  scale_color_manual(values = c("TRUE" = "#43006c", "FALSE" = "transparent"), guide = "none")
ggsave("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/TE_tree_long.pdf", height=23, width=9, device = cairo_pdf)
