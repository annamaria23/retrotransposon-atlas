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


my_cols <- paletteer::paletteer_d("ggthemes::Miller_Stone")
my_cols2 <- paletteer::paletteer_d("ggsci::planetexpress_futurama")

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

rm_sats_1 <- read_delim("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/TE_summary.txt") %>% 
    mutate(Species = sapply(strsplit(File, "/"), function(x) x[length(x) - 1])) %>% 
    left_join(full_sp)%>% 
    unique()%>% 
    na.omit()

rm_stats_dmel <- read_delim("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/repeatmasker_dmellib_summary.tsv", col_names=c("File", "Dmel")) %>% 
    mutate(Species = sapply(strsplit(File, "/"), function(x) x[length(x) - 2])) %>%  
    left_join(rm_sats_1, by = "Species") 

rm_sats <- read_delim("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/stats.csv") %>% 
    filter(pass=="PASS") %>% 
    rename(Species=species) %>% 
    left_join(rm_stats_dmel)

genome_lengths <- read_delim("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/genome_lengths.txt", col_names=c("Species", "GenomeSize")) %>% 
    left_join(full_sp)


#write.csv(rm_sats, "/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/rm_stats_combined.tsv", sep="\t")

sp_to_genus <- read_delim("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/Species_to_Genus.txt")

#Tree 
#tree <- read.iqtree("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/Kim_tree.tree")
tree <- read.iqtree("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/full_astral_fullAnnotation.tree")
labels <- tree@phylo$tip.label
labels <- tolower(labels)
labels <- gsub("_", ".", labels)
substr(labels, 1, 1) <- toupper(substr(labels, 1, 1))
tree@phylo$tip.label <- labels

tips_to_drop <- c(setdiff(tree@phylo$tip.label, rm_sats$SpeciesName), "D.pachea")
newtree <- drop.tip(tree@phylo, tips_to_drop)
p_tree <- ggtree(newtree, branch.length = "branch.length", layout = "fan", open.angle=15) + theme_tree()

Dmel_clade_species <- c("D.ficusphila", "D.eugracilis", "D.erecta", "D.orena", "D.teissieri", "D.santomea", "D.yakuba", "D.mauritiana",
                        "D.sechellia", "D.simulans", "D.mimetica", "D.takahashii", "D.pseudotakahashii", "D.biarmipes", "D.subpulchrella", "D.suzukii")

Zap_clade_species <- c("Z.bogoriensis", "Z.ornatus", "Z.vittiger", "Z.lachaisei", "Z.camerounensis", "Z.nigranus", "Z.capensis", "Z.davidi",
                      "Z.taronus", "Z.africanus", "Z.indianus", "Z.gabonicus", "Z.ghesquierei", "Z.inermis", "Z.kolodkinae", "Z.tuberculatus", "Z.tsacasi")

p_tree$data <- p_tree$data %>% mutate(dmel = label == "D.melanogaster")
p_tree$data <- p_tree$data %>% mutate(mel_clade = label %in% Dmel_clade_species)
p_tree$data <- p_tree$data %>% mutate(zap_clade = label %in% Zap_clade_species)

p_tree +
  #geom_hilight(node=getMRCA(as.phylo(p_tree), c(c("D.ficusphila", "D.santomea"))), fill="#F5F5DC", alpha=0.5, extendto=1)+
  #geom_hilight(node=getMRCA(as.phylo(p_tree), c(c("Z.bogoriensis", "Z.camerounensis"))), fill="#D7CE9FFF", alpha=0.5, extendto=1)+
  #geom_tiplab(size = 1, align = TRUE) + 
  #geom_tippoint(mapping = aes(subset = dmel), color = "red", size = 3, alpha = 1)+
  #geom_tippoint(mapping = aes(subset = mel_clade), color = "red", size = 3, alpha = 0.5)+
  geom_strip("D.ficusphila", "D.santomea", barsize=1.5, color="#638B66FF", offset=0.1, extend=0.5)+
  geom_strip("Z.bogoriensis", "Z.nigranus", barsize=1.5, color="#D7CE9FFF", offset=0.1, extend=0.5)+
  geom_tippoint(mapping = aes(subset = dmel), shape = 21, fill = "red", color = "black", size = 1, alpha = 1)+
  geom_fruit(
    data = sp_to_genus,
    geom = geom_col, 
    mapping = aes(y = SpeciesName, x = 1, fill = Genus),
    offset = 0.1, pwidth = 0.1,
    axis.params = list(axis="x", text.size = 1, text.angle = 0)
  )+
  scale_fill_manual(values = my_cols)+
  new_scale_fill()+
  geom_fruit(
    data = rm_sats %>% 
      mutate(Genome_Other = 100 - Retroelements - DNA_transposons) %>% 
      select(SpeciesName, Retroelements, DNA_transposons, Genome_Other) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent") %>%
      mutate(Category = ifelse(Category == "Genome_Other", NA, Category)),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 0.1, pwidth = 0.3,
    axis.params = list(axis="x", text.size = 1, text.angle = 0)
  ) +
  scale_fill_manual(values = c("#66A6FF", "#FF879F"), na.value = "grey80") +

  new_scale_fill() + 
  geom_fruit(
    data = rm_sats %>% 
      mutate(
        LINEs_pc = (LINEs / Retroelements+(1e-10)) * 100,
        LTR_pc = (LTR / Retroelements+(1e-10)) * 100,
        Retro_Other = 100 - LINEs_pc - LTR_pc
      ) %>%
      select(SpeciesName, LINEs_pc, LTR_pc, Retro_Other) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent") %>%
      mutate(Category = ifelse(Category == "Retro_Other", NA, Category)),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 0.1, pwidth = 0.3,
    axis.params = list(axis="x", text.size = 1, text.angle = 0)
  ) +
  scale_fill_manual(values = c("#AA4465", "#FFA69E"), na.value = "grey80")+

  new_scale_fill() + 
    geom_fruit(
    data = rm_sats %>% 
      mutate(
        Gypsy_pc = (Gypsy / LTR+(1e-10)) * 100,
        Pao_pc = (Pao / LTR+(1e-10)) * 100,
        Copia_pc = (Copia / LTR+(1e-10)) * 100,
        Other = 100 - Gypsy_pc - Pao_pc - Copia_pc
      ) %>%
      select(SpeciesName, Gypsy_pc, Pao_pc, Copia_pc, Other) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent") %>%
      mutate(Category = ifelse(Category == "Other", NA, Category)),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 0.1, pwidth = 0.3,
    axis.params = list(axis="x", text.size = 1, text.angle = 0)
  ) +
  scale_fill_manual(values = c("#DDA448", "#8CBCB9", "#BB342F"), na.value = "grey80") + 
  new_scale_fill() + 
    geom_fruit(
    data = rm_sats %>% 
      mutate(
        Dmel_PC = (Dmel / Gypsy+(1e-10)) * 100,
        Dmel_PC = ifelse(Dmel_PC>100, 100, Dmel_PC),
        Dmel_PC = ifelse(is.na(Gypsy), 0, Dmel_PC),
        Other = 100 - Dmel_PC
      ) %>%
      select(SpeciesName, Dmel_PC, Other) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent") %>%
      mutate(Category = ifelse(Category == "Other", NA, Category)),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 0.1, pwidth = 0.3,
    axis.params = list(axis="x", text.size = 1, text.angle = 0)
  ) +
scale_fill_manual(values = c("#777DA7"), na.value = "grey80")

ggsave("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/NEW_circular_tree_plot.pdf")


######### MELANOGASTER CLADE ONLY

subtree <- extract.clade(newtree, node = getMRCA(newtree, c("D.ficusphila", "D.santomea")))
p_tree_sub <- ggtree(subtree, branch.length = "branch.length", layout = "rectangular") + theme_tree2() +   geom_tiplab(aes(label = ""))

p_tree_sub +
  geom_tiplab(size = 3, align = TRUE) + 
  geom_fruit(
    data = rm_sats %>% 
      mutate(Genome_Other = 100 - Retroelements - DNA_transposons) %>% 
      select(SpeciesName, Retroelements, DNA_transposons, Genome_Other) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent") %>%
      mutate(Category = ifelse(Category == "Genome_Other", NA, Category)),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 1, pwidth = 1,
    axis.params = list(axis="x", text.size = 1, text.angle = -90)
  ) +
  scale_fill_manual(values = c("#66A6FF", "#FF879F"), na.value = "grey80") +
  new_scale_fill() + 
  geom_fruit(
    data = rm_sats %>% 
      mutate(
        LINEs_pc = (LINEs / Retroelements+(1e-10)) * 100,
        LTR_pc = (LTR / Retroelements+(1e-10)) * 100,
      ) %>%
      select(SpeciesName, LINEs_pc, LTR_pc) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent"),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 0.1, pwidth = 1,
    axis.params = list(axis="x", text.size = 3, text.angle = -90)
  ) +
  scale_fill_manual(values = c("#AA4465", "#FFA69E"), na.value = "grey80")+

  new_scale_fill() + 
    geom_fruit(
    data = rm_sats %>% 
      mutate(
        Gypsy_pc = (Gypsy / LTR+(1e-10)) * 100,
        Pao_pc = (Pao / LTR+(1e-10)) * 100,
        Copia_pc = (Copia / LTR+(1e-10)) * 100,
        Other = 100 - Gypsy_pc - Pao_pc - Copia_pc
      ) %>%
      select(SpeciesName, Gypsy_pc, Pao_pc, Copia_pc, Other) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent") %>%
      mutate(Category = ifelse(Category == "Other", NA, Category)),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 0.1, pwidth = 1,
    axis.params = list(axis="x", text.size = 1, text.angle = -90)
  ) +
  scale_fill_manual(values = c("#DDA448", "#8CBCB9", "#BB342F"), na.value = "grey80") + 
    new_scale_fill() + 
    geom_fruit(
    data = rm_sats %>% 
      mutate(
        Dmel_PC = (Dmel / Gypsy+(1e-10))*100,
        Dmel_PC = ifelse(Dmel_PC>100, 100, Dmel_PC),
        Dmel_PC = ifelse(is.na(Gypsy), 0, Dmel_PC),
        Other = 100 - Dmel_PC
      ) %>%
      select(SpeciesName, Dmel_PC, Other) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent") %>%
      mutate(Category = ifelse(Category == "Other", NA, Category)),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 0.1, pwidth = 1,
    axis.params = list(axis="x", text.size = 3, text.angle = -90)
  ) +
scale_fill_manual(values = c("#777DA7"), na.value = "grey80")


ggsave("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/mel_tree_plot.pdf", height=4, width=5)


subtree <- extract.clade(newtree, node = getMRCA(newtree, c("Z.bogoriensis", "Z.camerounensis")))
p_tree_sub <- ggtree(subtree, branch.length = "branch.length", layout = "rectangular") + theme_tree2() +   geom_tiplab(aes(label = ""))

p_tree_sub +
  geom_tiplab(size = 3, align = TRUE) + 
  geom_fruit(
    data = rm_sats %>% 
      mutate(Genome_Other = 100 - Retroelements - DNA_transposons) %>% 
      select(SpeciesName, Retroelements, DNA_transposons, Genome_Other) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent") %>%
      mutate(Category = ifelse(Category == "Genome_Other", NA, Category)),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 1, pwidth = 1,
    axis.params = list(axis="x", text.size = 1, text.angle = -90)
  ) +
  scale_fill_manual(values = c("#66A6FF", "#FF879F"), na.value = "grey80") +
  new_scale_fill() + 
  geom_fruit(
    data = rm_sats %>% 
      mutate(
        LINEs_pc = (LINEs / Retroelements+(1e-10)) * 100,
        LTR_pc = (LTR / Retroelements+(1e-10)) * 100,
      ) %>%
      select(SpeciesName, LINEs_pc, LTR_pc) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent"),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 0.1, pwidth = 1,
    axis.params = list(axis="x", text.size = 3, text.angle = -90)
  ) +
  scale_fill_manual(values = c("#AA4465", "#FFA69E"), na.value = "grey80")+

  new_scale_fill() + 
    geom_fruit(
    data = rm_sats %>% 
      mutate(
        Gypsy_pc = (Gypsy / LTR+(1e-10)) * 100,
        Pao_pc = (Pao / LTR+(1e-10)) * 100,
        Copia_pc = (Copia / LTR+(1e-10)) * 100,
        Other = 100 - Gypsy_pc - Pao_pc - Copia_pc
      ) %>%
      select(SpeciesName, Gypsy_pc, Pao_pc, Copia_pc, Other) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent") %>%
      mutate(Category = ifelse(Category == "Other", NA, Category)),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 0.1, pwidth = 1,
    axis.params = list(axis="x", text.size = 1, text.angle = -90)
  ) +
  scale_fill_manual(values = c("#DDA448", "#8CBCB9", "#BB342F"), na.value = "grey80") + 
    new_scale_fill() + 
    geom_fruit(
    data = rm_sats %>% 
      mutate(
        Dmel_PC = (Dmel / Gypsy+(1e-10))*100,
        Dmel_PC = ifelse(Dmel_PC>100, 100, Dmel_PC),
        Dmel_PC = ifelse(is.na(Gypsy), 0, Dmel_PC),
        Other = 100 - Dmel_PC
      ) %>%
      select(SpeciesName, Dmel_PC, Other) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent") %>%
      mutate(Category = ifelse(Category == "Other", NA, Category)),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 0.1, pwidth = 1,
    axis.params = list(axis="x", text.size = 3, text.angle = -90)
  ) +
scale_fill_manual(values = c("#777DA7"), na.value = "grey80")

ggsave("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/zap_tree_plot.pdf", height=4, width=5)


#Tree - full page for the supplementary

#tree <- read.iqtree("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/Kim_tree.tree")

tree <- read.iqtree("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/full_astral_fullAnnotation.tree")
labels <- tree@phylo$tip.label
labels <- tolower(labels)
labels <- gsub("_", ".", labels)
substr(labels, 1, 1) <- toupper(substr(labels, 1, 1))
tree@phylo$tip.label <- labels

tips_to_drop <- c(setdiff(tree@phylo$tip.label, rm_sats$SpeciesName), "D.pachea")
newtree <- drop.tip(tree@phylo, tips_to_drop)
p_tree <- ggtree(newtree, branch.length = "branch.length", layout = "rectangular") + theme_tree()

p_tree +
  geom_tiplab(size = 3, align = TRUE) + 
  geom_fruit(
    data = rm_sats %>% 
      mutate(Genome_Other = 100 - Retroelements - DNA_transposons) %>% 
      select(SpeciesName, Retroelements, DNA_transposons, Genome_Other) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent") %>%
      mutate(Category = ifelse(Category == "Genome_Other", NA, Category)),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 5, pwidth = 3,
    axis.params = list(axis="x", text.size = 2, text.angle = 0)
  ) +
  scale_fill_manual(values = c("#66A6FF", "#FF879F"), na.value = "grey80") +
  new_scale_fill() + 
  geom_fruit(
    data = rm_sats %>% 
      mutate(
        LINEs_pc = (LINEs / Retroelements+(1e-10)) * 100,
        LTR_pc = (LTR / Retroelements+(1e-10)) * 100,
        Retro_Other = 100 - LINEs_pc - LTR_pc
      ) %>%
      select(SpeciesName, LINEs_pc, LTR_pc, Retro_Other) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent") %>%
      mutate(Category = ifelse(Category == "Retro_Other", NA, Category)),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 0.5, pwidth = 3,
    axis.params = list(axis="x", text.size = 2, text.angle = 0)
  ) +
  scale_fill_manual(values = c("#AA4465", "#FFA69E"), na.value = "grey80")+

  new_scale_fill() + 
    geom_fruit(
    data = rm_sats %>% 
      mutate(
        Gypsy_pc = (Gypsy / LTR+(1e-10)) * 100,
        Pao_pc = (Pao / LTR+(1e-10)) * 100,
        Copia_pc = (Copia / LTR+(1e-10)) * 100,
        Other = 100 - Gypsy_pc - Pao_pc - Copia_pc
      ) %>%
      select(SpeciesName, Gypsy_pc, Pao_pc, Copia_pc, Other) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent") %>%
      mutate(Category = ifelse(Category == "Other", NA, Category)),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 0.5, pwidth = 3,
    axis.params = list(axis="x", text.size = 2, text.angle = 0)
  ) +
  scale_fill_manual(values = c("#DDA448", "#8CBCB9", "#BB342F"), na.value = "grey80") + 
  new_scale_fill() + 
    geom_fruit(
    data = rm_sats %>% 
      mutate(
        Dmel_PC = (Dmel / Gypsy+(1e-10)) * 100,
        Dmel_PC = ifelse(Dmel_PC>100, 100, Dmel_PC),
        Dmel_PC = ifelse(is.na(Gypsy), 0, Dmel_PC),
        Other = 100 - Dmel_PC
      ) %>%
      select(SpeciesName, Dmel_PC, Other) %>%
      pivot_longer(-SpeciesName, names_to = "Category", values_to = "Percent") %>%
      mutate(Category = ifelse(Category == "Other", NA, Category)),
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = Percent, fill = Category),
    offset = 0.5, pwidth = 3,
    axis.params = list(axis="x", text.size = 2, text.angle = 0)
  ) +
scale_fill_manual(values = c("#777DA7"), na.value = "grey80")+
  new_scale_fill() + 
  geom_fruit(
    data = genome_lengths,
    geom = geom_col,
    mapping = aes(y = SpeciesName, x = GenomeSize),
    offset = 0.5, pwidth = 5,
    axis.params = list(axis="x", text.size = 2, text.angle = 0),
    fill = "#618D48",
    alpha=0.7
  )
  ggsave("/mnt/scratchc/ghlab/annamaria/drosoTEs/blast_automation_attempt/FIGURES/NEW_LONG_tree_plot.pdf",height=20, width=8, device = cairo_pdf)
