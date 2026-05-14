# 02_assign_compartments_to_genes.R
# Map A/B compartments from bins to genes

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Load Hi-C compartments
hic <- read_tsv("data/derived/hic_compartments_250kb.tsv")

# Load gene annotation (from Salmobase)
gene <- read_tsv(
  "https://salmobase.org/datafiles/TSV/genes/AtlanticSalmon/Ssal_v3.1/Ensembl_genes.tsv"
)

# Example: assume chromosomes already match (chr1, chr2...)
gene_chr <- gene %>%
  rename(
    chrom = seqname,
    start = start,
    end = end
  ) %>%
  select(gene_id, chrom, start, end)

# Define bin size (must match earlier step!)
bin_size <- 250000

# Assign genes to bins using midpoint
gene_comp <- gene_chr %>%
  mutate(
    gene_mid = floor((start + end) / 2),
    bin_start = floor(gene_mid / bin_size) * bin_size,
    bin_end = bin_start + bin_size
  ) %>%
  left_join(
    hic,
    by = c("chrom" = "chrom", "bin_start" = "start", "bin_end" = "end")
  )

# Remove genes without assignment
gene_comp <- gene_comp %>%
  filter(!is.na(comp_gill), !is.na(comp_liver))

# Save result
write_tsv(gene_comp, "data/derived/gene_compartments.tsv")
