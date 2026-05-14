# 00_load_inputs.R
# Load libraries and input data

suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(GenomeInfoDb)
})

# Chromosome mapping
mapping <- read_tsv("data/raw/Ssal_v3.1_chromosomes.tsv")

# Gene annotation
gene <- read_tsv(
  "https://salmobase.org/datafiles/TSV/genes/AtlanticSalmon/Ssal_v3.1/Ensembl_genes.tsv"
)

# Repeat annotation
repeats <- read_tsv(
  "https://salmobase.org/datafiles/datasets/Salmon_repeats/Ssal_v3.1_repeats_5.2.bed",
  col_names = FALSE
)
colnames(repeats) <- c("chromosome", "start", "end", "annotation")

# Expression data
expr_data <- read.delim("data/derived/all_eve_genes.txt", header = TRUE)

