args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript plot_srna_length_rpm.R <input.bam>", call. = FALSE)
}

input_bam <- args[1]
output_pdf <- args[2]

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(tidyverse)
  library(patchwork)
})

message(paste("Processing: ", input_bam))

param <- ScanBamParam(what = c("strand"))
bam_data <- readGAlignments(input_bam, param = param)

# RPM = (Count / Total Aligned Reads) * 1,000,000
total_aligned_reads <- length(bam_data)
rpm_factor <- total_aligned_reads / 1e6

df_plot <- data.frame(
  length = width(bam_data),
  strand = as.character(strand(bam_data))
) %>%
  group_by(length, strand) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(rpm = count / rpm_factor) %>%
  mutate(rpm_adj = ifelse(strand == "-", -rpm, rpm))

final_plot <- ggplot(df_plot, aes(x = length, y = rpm_adj, fill = strand)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("+" = "#C492B1", "-" = "#3F7CAC")) +
  theme_classic(base_size=18) +
  scale_x_continuous(breaks = seq(min(df_plot$length), max(df_plot$length), 1)) +
  scale_y_continuous(labels = abs) +
  labs(
    x = "Read Length (nt)",
    y = "RPM",
    fill = "Strand"
  ) +
  theme(
  legend.position = "none",
  aspect.ratio = 1,
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
  axis.line = element_blank(),
  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)
)

ggsave(output_pdf, plot = final_plot, device = grDevices::cairo_pdf, width = 5, height = 5)

message(paste("Saved to:", output_pdf))
