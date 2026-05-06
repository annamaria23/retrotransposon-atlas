library(tidyverse)
library(ggplot2)
library(ggridges)
library(ggrepel)
library(forcats)


locs <- read_delim("~/hannonlab Dropbox/Anna-Maria Papameletiou/TE_annot/TE_insertion_locations/TEs_locations.txt", col_names = c("chr", "TE_start", "TE_end", "TE", "TE_score", "strand", "gene_chr", 
                                                                   "gene_start", "gene_end", "gene_id", "gene_score", "gene_strand",
                                                                   "distance")) %>% 
  mutate(Location = case_when(
    distance == -1 ~ "chr_without_genes",
    distance == 0  ~ "Intragenic",
    abs(distance) > 1 & abs(distance) <= 2000 ~ "Promoter (<=2000 bp)",
    TRUE ~ "Intergenic (>2000 bp)"
  )) %>% 
  mutate(
    Full_Length = ifelse(abs(TE_end - TE_start) > 4000, T, F)
  ) %>% 
  extract(TE, 
          into = c("Species", "TE"), 
          regex = "^(.*)_(.*)$", 
          remove = T)


locs %>%
  filter(Species == "Zafr") %>% 
  group_by(TE, Location, Full_Length) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) %>% 
  ggplot(data =., aes(y = TE, x = count, fill = Location)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  labs(
    x = "Full-length insertions",
    y = "TE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))+
  facet_wrap(~Full_Length, scales = "free")


locs %>% 
  filter(Species == "Dgun") %>% 
  ggplot(data = ., aes(x = TE_score, fill = Location))+
  geom_histogram(binwidth = 1, boundary = 0, color = "white") +
  theme_classic()

cols <- c("#247BA0", "#70C1B3", "#F25F5C", "#50514F")

locs %>% 
  filter(Species %in% c("Dyak", "Dsan", "Dtei", "Dore", "Dere", "Dsim", "Dsec", "Dmau", "Dmel",
                        "Dsuz", "Dsubp", "Dbia", "Dpse", "Dtak", "Dmim", "Deug", 
                        "Zafr", "Zbog", "Zcam", "Zcap", "Zdav", "Zgab", "Zghe", "Zind",
                        "Zine", "Zkol", "Zlac", "Znig", "Zorn", "Ztar", "Ztsa", "Ztub", "Zvit")) %>% 
  ggplot(data = ., aes(x = TE_score, fill = Location))+
  geom_histogram(binwidth = 2, boundary = 0, color = "white") +
  labs(x = "% divergence from consensus")+
  scale_fill_manual(values = cols)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        aspect.ratio = 1,
        text = element_text(size = 16))+
  facet_wrap(~Species, scales = "free_y", ncol=5)


locs %>% 
  filter(Species %in% c("Dyak", "Dsan", "Dtei", "Dore", "Dere", "Dsim", "Dsec", "Dmau", "Dmel",
                        "Dsuz", "Dsubp", "Dbia", "Dpse", "Dtak", "Dmim", "Deug", 
                        "Zafr", "Zbog", "Zcam", "Zcap", "Zdav", "Zgab", "Zghe", "Zind",
                        "Zine", "Zkol", "Zlac", "Znig", "Zorn", "Ztar", "Ztsa", "Ztub", "Zvit")) %>% 
  filter(Full_Length==T) %>% 
  ggplot(data = ., aes(x = TE_score, fill = Location))+
  geom_histogram(binwidth = 0.5, boundary = 0, color = "white") +
  labs(x = "% divergence from consensus")+
  scale_fill_manual(values = cols)+
  xlim(0,5)+
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        strip.background = element_blank(),
        strip.text = element_text(face = "italic"),
        aspect.ratio = 1,
        text = element_text(size = 16))+
  facet_wrap(~Species, scales = "free_y", ncol=5)


locs %>%
  filter(Species %in% c("Dyak", "Dsan", "Dtei", "Dore", "Dere", "Dsim", "Dsec", "Dmau", "Dmel",
                        "Dsuz", "Dsubp", "Dbia", "Dpse", "Dtak", "Dmim", "Deug", 
                        "Zafr", "Zbog", "Zcam", "Zcap", "Zdav", "Zgab", "Zghe", "Zind",
                        "Zine", "Zkol", "Zlac", "Znig", "Zorn", "Ztar", "Ztsa", "Ztub", "Zvit")) %>% 
  filter(TE %in% c("412", "blood", "mdg1", "Tabor", 
                   "Stalker", "Stalker4", "gypsy", "gtwin",
                   "gypsy2", "springer", "gypsy5", "Idefix",
                   "297", "Quasimodo","17.6", "Burdock", 
                   "micropia", "Circe", "chouto")) %>% 
  filter(Full_Length == T) %>% 
  group_by(TE, Location, Species) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count) * 100) %>% 
  ggplot(data =., aes(y = Species, x = count, fill = Location)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = cols)+
  theme_classic() +
  labs(
    x = "Full-length insertions",
    y = "TE") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text.y = element_text(angle = -30, vjust = 1, hjust = 1),
        strip.background = element_blank(),
        text = element_text(size = 15))+
  facet_wrap(~TE, scales = "free", ncol = 5)

locs %>%
  filter(Species %in% c("Dyak", "Dsan", "Dtei", "Dore", "Dere", "Dsim", "Dsec", "Dmau", "Dmel",
                        "Dsuz", "Dsubp", "Dbia", "Dpse", "Dtak", "Dmim", "Deug", 
                        "Zafr", "Zbog", "Zcam", "Zcap", "Zdav", "Zgab", "Zghe", "Zind",
                        "Zine", "Zkol", "Zlac", "Znig", "Zorn", "Ztar", "Ztsa", "Ztub", "Zvit")) %>% 
  filter(TE %in% c("412", "blood", "mdg1", "Tabor", 
                   "Stalker", "Stalker4", "gypsy", "gtwin",
                   "gypsy2", "springer", "gypsy5", "Idefix",
                   "297", "Quasimodo","17.6", "Burdock", 
                   "micropia", "Circe", "chouto")) %>% 
  filter(Full_Length == T) %>% 
  filter(TE_score <= 5) %>% 
  group_by(TE, Species) %>% 
  filter(n() >= 3) %>% 
  ungroup() %>%
  ggplot(data =., aes(y = Species, x = TE_score, fill = Species)) +
  geom_violin(alpha = 0.9, color = "black", scale = "width") +  
  geom_jitter(alpha=0.25, height = 0.3)+
  scale_fill_manual(values = my_colors) +
  theme_classic()+
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text.y = element_text(angle = -30, vjust = 1, hjust = 1),
        strip.background = element_blank(),
        text = element_text(size = 14),
        aspect.ratio = 1.25)+
  theme(legend.position = "none") +
  facet_wrap(~TE, scales = "free", ncol = 4)


locs %>%
  filter(Full_Length == TRUE) %>% 
  filter(TE_score < 1) %>% 
  group_by(TE, Species) %>%
  summarise(count = n(), .groups = "drop") %>%
  unique() %>% 
  mutate(TE = fct_reorder(TE, count, .fun = max)) %>% 
  ggplot(aes(y = TE, x = (count))) +
  geom_boxplot(outlier.shape = NA, color = "darkred", fill = "darkred", alpha=0.5, staplewidth = 0.8) + 
  geom_text_repel(
    data = . %>% filter(count > 25), 
    aes(label = Species),
    size = 4,
    box.padding = 0.2,
    point.padding = 0.1,
    direction = "x",
    segment.color = 'grey50'
  ) +
  geom_point(data = . %>% filter(count > 25), alpha = 0.5, size = 3, color = "darkred") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  labs(x = "Full-length insertions", y = "TE") +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.text.y = element_text(angle = -30, vjust = 1, hjust = 1),
        strip.background = element_blank(),
        text = element_text(size = 15),
        aspect.ratio = 1.5,
        panel.grid.major.y = element_line(colour = "lightgrey", linewidth = 0.5),
        panel.grid.minor.y = element_line(colour = "lightgrey", linewidth = 0.25))

  
