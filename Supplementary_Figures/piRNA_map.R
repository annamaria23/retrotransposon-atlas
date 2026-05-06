#!/usr/bin/env Rscript

#conda activate srna_plots

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: No BAM file provided. Usage: Rscript plot_mirrored_logo.R <input.bam>", call. = FALSE)
}

input_bam <- args[1]
output_pdf <- args[2]

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(ggseqlogo)
  library(tidyverse)
  library(Biostrings)
  library(patchwork)
})

message(paste("Processing:", input_bam))

param <- ScanBamParam(what = c("seq", "strand"))
bam_data <- readGAlignments(input_bam, param = param)

raw_seqs <- mcols(bam_data)$seq
strands <- as.character(strand(bam_data))

logo_len <- 15
corrected_seqs <- raw_seqs
minus_idx <- which(strands == "-")
corrected_seqs[minus_idx] <- reverseComplement(raw_seqs[minus_idx])

create_rna_matrix <- function(seq_set_dna, length_cutoff) {
  valid_seqs <- seq_set_dna[width(seq_set_dna) >= length_cutoff]
  
  if(length(valid_seqs) == 0) return(NULL)
  valid_seqs %>%
    subseq(start=1, end=length_cutoff) %>%
    as.character() %>%
    stringr::str_replace_all("T", "U") %>%
    RNAStringSet() %>%
    consensusMatrix(as.prob = TRUE) %>%
    .[c("A","C","G","U"), ]
}

sense_mat <- create_rna_matrix(corrected_seqs[strands == "+"], logo_len)
anti_mat  <- create_rna_matrix(corrected_seqs[strands == "-"], logo_len)

p <- ggseqlogo(anti_mat, method = 'bits') +
  #scale_y_reverse() +
  theme_classic(base_size = 18) +
  theme(plot.margin = margin(t = 5)) + 
  labs(y = "Bits", x = "Position (nt)")+
  theme(
    aspect.ratio = 0.5,
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 24),
    axis.text = element_text(color = "black", size = 16),
    axis.title = element_text(size = 20),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave(output_pdf, 
       plot = p, 
       device = grDevices::cairo_pdf, 
       width = 6, height = 3.5)

message(paste("Saved to:", output_pdf))
