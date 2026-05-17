library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

rm_full_stats <- read_delim("~/hannonlab Dropbox/Anna-Maria Papameletiou/TE_annot/Figure1_General/rm_stats_combined.tsv")

len %>% 
  pivot_longer(
    cols = c(GenomeLength, Total.len), 
    names_to = "Len_Type", 
    values_to = "Length"
  ) %>% 
ggplot(data = ., aes(y = Length, x= sp, fill = Len_Type))+
  geom_violin()+
  geom_jitter(alpha=0.5, color = "grey80")+
  theme_classic()+
  scale_color_manual(values = (my_cols))+
    theme(
    #aspect.ratio = 1,
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(color = "black")
  )


rm_full_stats %>% 
  pivot_longer(
    cols = c(Total, Retroelements, DNA_transposons), 
    names_to = "Len_Type", 
    values_to = "Length"
  ) %>% 
  mutate(Len_Type = factor(Len_Type, levels = c("DNA_transposons", "Total", "Retroelements"))) %>%
ggplot(data = ., aes(y = Length, x= Len_Type, fill = Len_Type, , color = Len_Type))+
  geom_line(aes(group = SpeciesName), color = "grey50", alpha = 0.2) + 
  geom_boxplot(alpha=0.5, outlier.shape = NA, staplewidth = 1)+
  scale_color_manual(values = c("#67a6ff", "grey40", "#ff879f"))+
  scale_fill_manual(values = c("#67a6ff", "grey40", "#ff879f"))+
  geom_point(alpha=0.1, color = "black")+
  labs(x = "TE", y = "% Genome")+
  theme_classic(base_size = 16)+
    theme(
    aspect.ratio = 1,
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(color = "black")
  )

len %>% 
  pivot_longer(
    cols = c(LTR, LINEs), 
    names_to = "Len_Type", 
    values_to = "Length"
  ) %>% 
ggplot(data = ., aes(y = Length, x= Len_Type, fill = Len_Type, , color = Len_Type))+
  geom_line(aes(group = SpeciesName), color = "grey50", alpha = 0.2) + 
  geom_boxplot(alpha=0.5, outlier.shape = NA, staplewidth = 1)+
  scale_color_manual(values = c("#ffa69e", "#aa4465"))+
  scale_fill_manual(values = c("#ffa69e", "#aa4465"))+
  geom_point(alpha=0.1, color = "black")+
  labs(x = "TE", y = "% Genome")+
  theme_classic(base_size = 16)+
    theme(
    aspect.ratio = 1,
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(color = "black")
  )


len %>% 
  pivot_longer(
    cols = c(Copia, Gypsy, Pao), 
    names_to = "Len_Type", 
    values_to = "Length"
  ) %>% 
ggplot(data = ., aes(y = Length, x= Len_Type, fill = Len_Type, , color = Len_Type))+
  geom_line(aes(group = SpeciesName), color = "grey50", alpha = 0.2) + 
  geom_boxplot(alpha=0.5, outlier.shape = NA, staplewidth = 1)+
  scale_color_manual(values = c("#dda448", "#8cbcb9", "#bb342f"))+
  scale_fill_manual(values = c("#dda448", "#8cbcb9", "#bb342f"))+
  geom_point(alpha=0.1, color = "black")+
  labs(x = "TE", y = "% Genome")+
  theme_classic(base_size = 16)+
    theme(
    aspect.ratio = 1,
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(color = "black")
  )


rm_full_stats %>% 
  mutate(Gypsy_Dmel = (ifelse(Dmel>Gypsy, Gypsy, Dmel))) %>% 
  pivot_longer(
    cols = c(Gypsy, Gypsy_Dmel), 
    names_to = "Len_Type", 
    values_to = "Length"
  ) %>% 
ggplot(data = ., aes(y = Length, x= Len_Type, fill = Len_Type, , color = Len_Type))+
  geom_line(aes(group = SpeciesName), color = "grey50", alpha = 0.2) + 
  geom_boxplot(alpha=0.5, outlier.shape = NA, staplewidth = 1)+
  scale_color_manual(values = c("#8cbcb9", "#777DA7"))+
  scale_fill_manual(values = c("#8cbcb9", "#777DA7"))+
  geom_point(alpha=0.1, color = "black")+
  labs(x = "TE", y = "% Genome")+
  theme_classic(base_size = 16)+
    theme(
    aspect.ratio = 1,
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(color = "black")
  )
