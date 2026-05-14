# 01_assign_compartments_to_bins.R
# Load Hi-C eigenvectors and define A/B compartments

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Load eigenvector files
HiCG <- read_tsv("data/derived/HiC_G.250000_compartments.cis.vecs.tsv")
HiCL <- read_tsv("data/derived/HiC_L.250000_compartments.cis.vecs.tsv")

# Rename columns for clarity
HiCG <- HiCG %>%
  rename(
    E1_gill = E1
  )

HiCL <- HiCL %>%
  rename(
    E1_liver = E1
  )

# Merge and assign compartments
hic <- HiCG %>%
  left_join(HiCL, by = c("chrom", "start", "end")) %>%
  mutate(
    comp_gill  = ifelse(E1_gill  > 0, "A", "B"),
    comp_liver = ifelse(E1_liver > 0, "A", "B")
  )

# Save clean output
write_tsv(hic, "data/derived/hic_compartments_250kb.tsv")
