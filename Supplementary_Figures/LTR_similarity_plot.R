library(tidyverse)
library(readr)
library(ggplot2)
library(patchwork)

ltr_sim <- read_delim("LTR_similarity.txt", delim = " ", col_names = c("LTR_1", "LTR_2", "Sim")) %>% 
  unique() %>% 
  filter(LTR_1 != LTR_2) %>% 
  mutate(Sim = as.numeric(str_remove(Sim, "%"))) %>% 
  separate_wider_regex(
    LTR_1, 
    patterns = c(
      Name = ".*?",
      "::.*?:",
      Start = "[0-9]+", 
      "-",
      End = "[0-9]+" 
    )
  ) %>%
  mutate(across(c(Start, End), as.numeric)) %>% 
  separate_wider_regex(
    Name, 
    patterns = c(
      Species = ".*", 
      "_", 
      TE = ".*" 
    )
  ) %>% 
  filter(Start == 0) %>% 
  filter(End >=150)

ltr_summary <- ltr_sim %>%
  group_by(TE) %>%
  summarize(
    pct_high = mean(Sim > 95) * 100,
    median_sim = median(Sim)
  ) %>%
  mutate(TE = fct_reorder(TE, median_sim)) %>% 
  mutate(TE = fct_reorder(TE, median_sim))

p1 <- ggplot(ltr_sim, aes(x = fct_reorder(TE, Sim, .fun = median), y = Sim)) +
  geom_violin(fill = "steelblue", color = "steelblue", alpha = 0.5) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  coord_flip(ylim = c(80, 100)) + 
  labs(x = "TE Family", y = "Percent Similarity") +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank() 
  )

p2 <-ggplot(ltr_summary, aes(x = pct_high, y = TE)) +
  geom_col(aes(x = 100), fill = "grey90") +
  geom_col(aes(fill = pct_high), fill = "steelblue") +
  geom_text(aes(label = paste0(round(pct_high), "%")), 
            hjust = 1.1, color = "white", size = 4, fontface = "bold") +
  labs(x = "% LTRs > 95 % similarity", y = NULL) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.y = element_blank(),
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    axis.line = element_blank(),
    axis.ticks.y = element_blank()
  )

p1 + p2 + plot_layout(widths = c(2, 1))

ggsave(
  filename = "LTR_similarities.pdf", 
  device = cairo_pdf,
  width = 6,
  height = 14)

