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
library(stringr)

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

full_sp <- rbind(species, more_species)

# Data 
SW_scores <- read_delim("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/ALL_sw_scores.tsv", 
col_names = c("Species", "TF", "Similarity", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "SW")) %>%  
select(c("Species", "TF", "Similarity", "SW")) %>%  
mutate(Species = sub("_.*", "", Species)) %>%
group_by(Species, TF) %>%
filter(SW == max(SW)) %>%
distinct(Species, TF, .keep_all = TRUE) %>%
ungroup()%>%
filter(Similarity >= 50) %>%
left_join(full_sp) %>%
complete(SpeciesName, TF) 

busco_scores <- read_delim("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/busco_summary.tsv", 
col_names=c("Score", "Species")) %>%
mutate(Species = str_split_i(Species, "/", 9))%>% 
left_join(full_sp) %>%
unique() %>%
mutate(score_bin = case_when(Score < 85 ~ "<85",
  Score >=85 & Score <=95 ~ "85-95",
  Score > 95 ~ ">95"))

tj_sw <- SW_scores %>% filter(TF == "tj") %>% mutate(SW = SW/555)
CHES_sw <- SW_scores %>% filter(TF == "CHES-1-like") %>% mutate(SW = SW/1268)
CrebB_sw <- SW_scores %>% filter(TF == "CrebB") %>% mutate(SW = SW/331)
Cf2_sw <- SW_scores %>% filter(TF == "Cf2") %>% mutate(SW = SW/537)
GATAd_sw <- SW_scores %>% filter(TF == "GATAd") %>% mutate(SW = SW/842)
BEAF_sw <- SW_scores %>% filter(TF == "BEAF-32") %>% mutate(SW = SW/283)

somatic_sw <- rbind(tj_sw, CHES_sw)
germline_sw <- rbind(CrebB_sw, Cf2_sw, GATAd_sw, BEAF_sw)


# Tree
#tree <- read.iqtree("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/Kim_tree.tree")
tree <- read.iqtree("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/full_astral_fullAnnotation.tree")
labels <- tree@phylo$tip.label
labels <- tolower(labels)
labels <- gsub("_", ".", labels)
substr(labels, 1, 1) <- toupper(substr(labels, 1, 1))
tree@phylo$tip.label <- labels

tips_to_drop <- c(setdiff(tree@phylo$tip.label, SW_scores$SpeciesName), "D.pachea")
newtree <- drop.tip(tree@phylo, tips_to_drop)
p_tree <- ggtree(newtree, branch.length = "branch.length", layout = "fan", open.angle=90) + theme_tree()
p_tree <- rotate_tree(p_tree, 90)

# Data
'''
p_tree +
  geom_fruit(
    data = busco_scores,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = "", fill = score_bin),
    offset = 0.1, pwidth = 4,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1))+
  scale_fill_manual(values = c(">95"="#dadad6", "85-95"="#9f9f9d", "<85"="#4A4A48"), name = "BUSCO completeness %")+
  new_scale_fill() + 
  geom_fruit(
    data = tj_sw,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TF, fill = as.numeric(SW)),
    offset = 0.1, pwidth = 4,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1))+
  scale_fill_gradient(low = "#bed0ff", high = "#00299a", name = "tj SW score", na.value="white") + 
  new_scale_fill() + 
  geom_fruit(
    data = CHES_sw,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TF, fill = as.numeric(SW)),
    offset = 0.1, pwidth = 4,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1))+
  scale_fill_gradient(low = "#bed0ff", high = "#00299a", name = "CHES-1-like SW score", na.value="white") + 
  new_scale_fill() + 
    geom_fruit(
    data = CrebB_sw,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TF, fill = as.numeric(SW)),
    offset = 0.1, pwidth = 4,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1)) +
    scale_fill_gradient(low = "#ffb5db", high = "#6a0137", name = "CrebB SW score", na.value="white")+
  new_scale_fill() + 
    geom_fruit(
    data = Cf2_sw,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TF, fill = as.numeric(SW)),
    offset = 0.1, pwidth = 4,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1))+
    scale_fill_gradient(low = "#ffb5db", high = "#6a0137", name = "Cf2 SW score", na.value="white")+
  new_scale_fill() + 
    geom_fruit(
    data = GATAd_sw,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TF, fill = as.numeric(SW)),
    offset = 0.1, pwidth = 4,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1))+
    scale_fill_gradient(low = "#ffb5db", high = "#6a0137", name = "GATAd TFs", na.value="white")+
  new_scale_fill() + 
    geom_fruit(
    data = BEAF_sw,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TF, fill = as.numeric(SW)),
    offset = 0.1, pwidth = 4,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1))+
    scale_fill_gradient(low = "#ffb5db", high = "#6a0137", name = "BEAF-32 TFs", na.value="white")
'''

p_tree$data <- p_tree$data %>% mutate(dmel = label == "D.melanogaster")

p_tree +
  geom_tippoint(mapping = aes(subset = dmel), shape = 21, fill = "red", color = "black", size = 3, alpha = 1)+
  geom_fruit(
    data = busco_scores,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = "", fill = score_bin),
    offset = 0.1, width = 4,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1))+
  scale_fill_manual(values = c(">95"="#dadad6", "85-95"="#9f9f9d", "<85"="#4A4A48"), name = "BUSCO completeness %")+
  new_scale_fill() + 
  geom_fruit(
    data = somatic_sw,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TF, fill = as.numeric(SW)),
    offset = 0.1, width = 4,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1))+
  scale_fill_gradient(low = "#bed0ff", high = "#00299a", name = "tj SW score", na.value="white") + 
  new_scale_fill() + 
  geom_fruit(
    data = germline_sw,
    geom = geom_tile,
    mapping = aes(y = SpeciesName, x = TF, fill = as.numeric(SW)),
    offset = 0.1, pwidth = 0.6, width=4,
    axis.params = list(axis = "x", text.size = 3 , text.angle = 0, hjust=-0.1)) +
    scale_fill_gradient(low = "#ffb5db", high = "#6a0137", name = "CrebB SW score", na.value="white")
ggsave("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/TF_SW.pdf", device = cairo_pdf, height=8, width=8)



