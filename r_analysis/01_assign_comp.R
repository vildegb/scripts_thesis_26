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

# 30% thresholds
q_low_gill   <- quantile(hic$E1_gill,  0.30, na.rm = TRUE)
q_high_gill  <- quantile(hic$E1_gill,  0.70, na.rm = TRUE)

q_low_liver  <- quantile(hic$E1_liver, 0.30, na.rm = TRUE)
q_high_liver <- quantile(hic$E1_liver, 0.70, na.rm = TRUE)

hic <- hic %>%
  mutate(
    comp30_gill = case_when(
      E1_gill <= q_low_gill  ~ "B",
      E1_gill >= q_high_gill ~ "A",
      TRUE ~ "Ambiguous"
    ),
    comp30_liver = case_when(
      E1_liver <= q_low_liver  ~ "B",
      E1_liver >= q_high_liver ~ "A",
      TRUE ~ "Ambiguous"
    )
  )

# Save clean output
write_tsv(hic, "data/derived/hic_compartments_250kb.tsv")
