# 05_upset_compartments.R
# UpSet plot for duplicated genes across A/B compartments

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(UpSetR)
})

# Load data

final_table <- read_csv("data/derived/final_table.csv")

# Keep one row per gene
final_genes <- final_table %>%
  distinct(Ssal, .keep_all = TRUE)

# Prepare duplicated gene set

gene_comp <- final_genes %>%
  filter(dup_type == "duplicate") %>%
  filter(!is.na(comp_gill), !is.na(comp_liver)) %>%
  group_by(Ssal) %>%
  summarise(
    comp_gill  = first(comp_gill),
    comp_liver = first(comp_liver),
    .groups = "drop"
  )

# Create presence/absence matrix

gene_comp <- gene_comp %>%
  mutate(
    A_gill  = comp_gill  == "A",
    B_gill  = comp_gill  == "B",
    A_liver = comp_liver == "A",
    B_liver = comp_liver == "B"
  )

upset_matrix <- data.frame(
  A_gill  = as.numeric(gene_comp$A_gill),
  B_gill  = as.numeric(gene_comp$B_gill),
  A_liver = as.numeric(gene_comp$A_liver),
  B_liver = as.numeric(gene_comp$B_liver)
)


# Plot

png(
  filename = "figures/Fig_upset_compartments.png",
  width = 2000,
  height = 1500,
  res = 300
)

upset(
  upset_matrix,
  sets = c("A_gill", "B_gill", "A_liver", "B_liver"),
  order.by = "freq",
  mainbar.y.label = "Number of duplicated genes",
  sets.x.label = "Compartment-Tissue combinations",
  text.scale = c(1.2, 1.2, 1, 1, 1, 1.6)
)
