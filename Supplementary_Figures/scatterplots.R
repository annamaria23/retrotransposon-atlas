library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

my_cols <- paletteer_d("ggthemes::Miller_Stone")
my_cols2 <- paletteer_d("ggsci::planetexpress_futurama")

len <- read_delim("~/hannonlab Dropbox/Anna-Maria Papameletiou/TE_annot/Figure1_General/genome_lengths.txt", col_names = c("Species", "GenomeLength")) %>% 
  left_join(rm_sats) %>% 
  na.omit() %>% 
  mutate(Total.len = (Total*GenomeLength)/100,
         Retroelements.len = (Retroelements*GenomeLength)/100,
         DNA_transposons.len = (DNA_transposons*GenomeLength)/100,
         LTR.len = (LTR*GenomeLength)/100,
         ) %>% 
  mutate(sp = substr(SpeciesName, 1, 1)) %>% 
  mutate(genus = case_when(
    sp == "A" ~ "Aedes",
    sp == "C" ~ "Chymomyza",
    sp == "D" ~ "Drosophila",
    sp == "H" ~ "Hirtodrosophila",
    sp == "L" ~ "Lordiphosa",
    sp == "S" ~ "Scaptodrosophila",
    sp == "Z" ~ "Zaprionus"
  )) %>% 
  unique()

len <- read_delim("~/hannonlab Dropbox/Anna-Maria Papameletiou/TE_annot/Figure1_General/genome_lengths.txt", col_names = c("Species", "GenomeLength")) %>% 
  left_join(rm_sats) %>% 
  na.omit() %>% 
  mutate(Total.len = (Total*GenomeLength)/100,
         Retroelements.len = (Retroelements*GenomeLength)/100,
         DNA_transposons.len = (DNA_transposons*GenomeLength)/100,
         LTR.len = (LTR*GenomeLength)/100,
         ) %>% 
  mutate(sp = substr(SpeciesName, 1, 1)) %>% 
  left_join(genus) %>% 
  mutate(Genus = ifelse(is.na(Genus), "Drosophila", Genus)) %>% 
  unique()

sp_genus <- len %>% 
  dplyr::select(SpeciesName, Genus)


# Plot all v all
######

ggplot(data = len, aes(x = GenomeLength, y = Total, color = Genus))+
  geom_point(alpha=0.9)+
  geom_text_repel(data = subset(len, Genus %in% c("Chymomyza", "Amiota", "Cacoxenus")), 
                  aes(label = SpeciesName),
                  size = 5, 
                  nudge_x = 20000000
                  )+
    geom_text_repel(data = subset(len, SpeciesName %in% c("D.melanogaster")), 
                  aes(label = SpeciesName),
                  size = 5, 
                  #nudge_x = 100000000,
                  nudge_y = 5
                  )+
  theme_classic(base_size = 16)+
  stat_cor(aes(group = 1), method = "pearson")+
  scale_color_manual(values = (my_cols))+
  theme(
    aspect.ratio = 1,
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(color = "black")
  )+
  labs(y = "Total Transposons (%)",
      x = "Total Genome Size (bp)")

ggplot(data = len, aes(x = GenomeLength, y = DNA_transposons, color = Genus))+
  geom_point(alpha=0.9)+
  geom_text_repel(data = subset(len, Genus %in% c("Chymomyza", "Amiota", "Cacoxenus")), 
                  aes(label = SpeciesName),
                  size = 5, 
                  nudge_x = 20000000)+
    geom_text_repel(data = subset(len, SpeciesName %in% c("D.melanogaster")), 
                  aes(label = SpeciesName),
                  size = 5, 
                  nudge_y = -5)+
  theme_classic(base_size = 16)+
  stat_cor(aes(group = 1), method = "pearson")+
  scale_color_manual(values = (my_cols))+
  theme(
    aspect.ratio = 1,
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(color = "black")
  )+
  labs(y = "DNA Transposons (%)",
      x = "Total Genome Size (bp)")

ggplot(data = len, aes(x = GenomeLength, y = Retroelements, color = Genus))+
  geom_point(alpha=0.9)+
  geom_text_repel(data = subset(len, Genus %in% c("Chymomyza", "Amiota", "Cacoxenus")), 
                  aes(label = SpeciesName),
                  size = 5, 
                  nudge_x = 20000000)+
    geom_text_repel(data = subset(len, SpeciesName %in% c("D.melanogaster")), 
                  aes(label = SpeciesName),
                  size = 5, 
                  nudge_y = -5)+
  theme_classic(base_size = 16)+
  stat_cor(aes(group = 1), method = "pearson")+
  scale_color_manual(values = (my_cols))+
  theme(
    aspect.ratio = 1,
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(color = "black")
  )+
  labs(y = "RNA Transposons (%)",
      x = "Total Genome Size (bp)")


ggplot(data = len, aes(x = GenomeLength, y = LTR, color = Genus))+
  geom_point(alpha=0.9)+
  geom_text_repel(data = subset(len, Genus %in% c("Chymomyza", "Amiota", "Cacoxenus")), 
                  aes(label = SpeciesName),
                  size = 5, 
                  nudge_x = 20000000)+
    geom_text_repel(data = subset(len, SpeciesName %in% c("D.melanogaster")), 
                  aes(label = SpeciesName),
                  size = 5, 
                  nudge_y = 5)+
  theme_classic(base_size = 16)+
  stat_cor(aes(group = 1), method = "pearson")+
  scale_color_manual(values = (my_cols))+
  theme(
    aspect.ratio = 1,
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(color = "black")
  )+
  labs(y = "LTR Retrotransposons (%)",
      x = "Total Genome Size (bp)")
