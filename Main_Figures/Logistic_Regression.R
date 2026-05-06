library(tidyverse)
library(ggplot2)
library(dplyr)
library(broom)
library(broom.mixed)
library(lme4)
library(reshape2)
library(proxy)

sc_expr <- read_delim("FlyCellAtlas_slimmed_gene_expression_fb_2025_04.fmt.tsv") # FCA gene expression, downloaded from FlyBase
fly_tfs <- read_delim("flymine_results_2025-10-22T13-50-21.tsv", col_names = c("geneID", "gene_Symbol", "geneFbID", "Dmel")) # TFs in D. melanogaster, downloaded from FlyMine

tf_expr_all <- left_join(fly_tfs, sc_expr) %>% 
  separate(female_germline_cell, into = c("female_germline_cell_expr", "female_germline_cell_cells"), sep = ":", convert = TRUE) %>% 
  mutate(female_germline_cell_cells = gsub("%", "", female_germline_cell_cells)) %>% 
  separate(somatic_cell_of_ovariole, into = c("somatic_cell_of_ovariole_expr", "somatic_cell_of_ovariole_cells"), sep = ":", convert = TRUE) %>% 
  mutate(somatic_cell_of_ovariole_cells = gsub("%", "", somatic_cell_of_ovariole_cells)) %>% 
  mutate(somatic_cell_of_ovariole_cells = as.numeric(somatic_cell_of_ovariole_cells),
         female_germline_cell_cells = as.numeric(female_germline_cell_cells)) %>% 
  mutate(tf_type = case_when(
    female_germline_cell_cells >= 20 & somatic_cell_of_ovariole_cells >= 20 ~ "both",
    female_germline_cell_cells >= 20 ~ "germline",
    somatic_cell_of_ovariole_cells >= 20 ~ "somatic",
    TRUE ~ "none"))

gl_tfs <- unique(tf_expr_all %>% filter(tf_type == "germline")) %>% pull(gene_Symbol)
soma_tfs <- unique(tf_expr_all %>% filter(tf_type == "somatic")) %>% pull(gene_Symbol) 
both_tfs <- unique(tf_expr_all %>% filter(tf_type == "both")) %>% pull(gene_Symbol)
ovary_tfs <- c(gl_tfs, soma_tfs, both_tfs) %>% unique()

env_status <- read_delim("ENV_locs.bed", col_names = c("sequence_name", "start", "end", "ENV", "p", "strand")) %>% 
  filter(end - start > 900) %>% 
  filter(end - start < 2000) 

fimo_output <- read_delim("all_TFs_LTRs_fimo.tsv") %>% 
  mutate(motif_alt_id = word(motif_alt_id, start = 3, end = -1, sep = fixed("."))) %>% 
  select(c(motif_alt_id, sequence_name)) %>% 
  filter(motif_alt_id %in% ovary_tfs) %>% 
  mutate(motif_alt_id = as.factor(motif_alt_id)) %>% 
  unique() %>% 
  group_by(motif_alt_id) %>%
  filter(n() >= 30) %>%
  ungroup() %>% 
  left_join(env_status, relationship = "many-to-many") %>% 
  mutate(y = ifelse(is.na(ENV), F, T)) %>% 
  select(motif_alt_id, sequence_name, y) %>% 
  filter(!str_detect(sequence_name, "gypsy12")) %>% 
  filter(!str_detect(sequence_name, "Bica")) %>% 
  filter(!str_detect(sequence_name, "Transpac")) %>% 
  filter(!str_detect(sequence_name, "McClintock")) %>% 
  filter(!str_detect(sequence_name, "gypsy9")) %>%  # Removing TEs with few species or too few ENV hits
  unique() 

df_wide_motifs <- fimo_output %>%
  mutate(motif_present = 1) %>%
  pivot_wider(
    id_cols = sequence_name,
    names_from = motif_alt_id,
    values_from = motif_present,
    values_fn = max, values_fill = 0)

data_for_glmm <- df_wide_motifs %>% 
  mutate(te_id = factor(str_extract(sequence_name, "[^_]+$"))) %>%
  #mutate(te_id = factor(str_remove(sequence_name, "^[_^]+")))  %>% 
  left_join(df_outcome, by = "sequence_name") %>%
  select(-sequence_name)

model <- glmer(y_outcome ~ . - te_id + (1 | te_id), data = data_for_glmm, family = binomial(link = "logit"),
               control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

# summary(model)


#### Plot the model outputs ####

df <- tidy(model) %>%
  filter(term != "(Intercept)") %>% 
  filter(term != "sd__(Intercept)") %>% 
  mutate(
    p_adj = p.adjust(p.value, method = "fdr"),
    OR = exp(estimate),                  # odds ratio
    sig = p_adj < 0.05,                # significance flag
    term = gsub("motif_alt_id", "", term),
    term = gsub("`", "", term),
    term = reorder(term, OR)
  ) %>% 
  mutate(type = ifelse(term %in% soma_tfs, "soma", NA),
         type = ifelse(term %in% gl_tfs, "germline", type),
         type = ifelse(term %in% both_tfs, "both", type)) %>% 
  mutate(sig = ifelse(sig == TRUE, type, sig))

ggplot(df, aes(x = term, y = OR, color = sig)) +
  geom_segment(aes(x = term, xend = term, y = 1, yend = OR), 
               color = "gray90", size = 1) +
  geom_errorbar(aes(ymin = exp(estimate - 1.96 * std.error),
                    ymax = exp(estimate + 1.96 * std.error)),
                width = 0.3, size = 1.2, alpha = 0.8) +
  geom_point(aes(size = -log10(p.value)), alpha = 1) +
  scale_y_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10), 
                labels = c("0.1", "0.5", "1", "2", "5", "10")) +
  coord_flip() +
  scale_color_manual(values = c("soma" = "#AEC3FF", 
                                "both" = "#AE9CEC", 
                                "germline" = "#FF66B6", 
                                "FALSE" = "gray80"),
                     labels = c("FALSE" = "Not significant", 
                                "soma" = "Soma", 
                                "germline" = "Germline", 
                                "both" = "Both"),
                     name = "Tissue Type") +
  scale_size_continuous(range = c(3, 8), name = "-log10(p)") +
  labs(x = NULL, y = "Odds Ratio (95% CI)") +
  geom_hline(yintercept = 1, lty = "dashed", color = "gray40", size = 0.8) +
  theme_minimal(base_size = 20) + 
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(size = 22),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.margin = margin(10, 20, 10, 10)
)

#### Example Venn Diagram ####

venn_data <- fimo_output %>% 
  mutate(TE = factor(str_extract(sequence_name, "[^_]+$"))) %>%
  mutate(Species = factor(str_remove(sequence_name, "^[_^]+"))) %>% 
  group_by(sequence_name) %>%
  summarise(
    CrebB = "CrebB" %in% motif_alt_id,
    ENV = first(y)
  ) %>% 
  select(-sequence_name)

fit <- euler(venn_data)

pdf("CrebB_venn.pdf", height = 2, width=2, useDingbats = F)
plot(fit,
     quantities = TRUE, 
     fills = c("#AE9CEC", "white"), 
     alpha = 1,
     edges = "black",
)
dev.off()


#### Plot the Jaccard Index Heatmap ### 

binary_matrix <- fimo_output %>%
  mutate(present = 1) %>%
  filter(motif_alt_id %in% c("tj", "CHES-1-like", "Cf2", "BEAF-32", "GATAd", "CrebB")) %>% 
  pivot_wider(id_cols = sequence_name, 
              names_from = motif_alt_id, 
              values_from = present, 
              values_fill = 0) %>%
  column_to_rownames("sequence_name") %>%
  as.matrix()

j_dist <- proxy::dist(t(binary_matrix), method = "Jaccard")
j_matrix <- 1 - as.matrix(j_dist)

hc <- hclust(as.dist(1 - j_matrix))
motif_order <- hc$labels[hc$order]
j_matrix <- j_matrix[motif_order, motif_order]

j_matrix[upper.tri(j_matrix, diag = TRUE)] <- NA

plot_df <- melt(j_matrix, na.rm = TRUE) %>%
  mutate(Var1 = factor(Var1, levels = motif_order),
         Var2 = factor(Var2, levels = motif_order))

ggplot(plot_df, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "white", high = "grey20", 
                      name = "Jaccard Index", na.value = "white") +
  coord_fixed() +
  labs(title = "Motif Co-occurrence", x = NULL, y = NULL) +
  theme_minimal(base_size = 22) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black", face = "bold"),
    plot.title = element_text(face = "bold", size = 26, margin = margin(b=15)),
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 16),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20)
  )









