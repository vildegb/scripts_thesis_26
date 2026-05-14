# 04_compartment_switching.R
# Compartment switching between gill and liver

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# Load data
final_table <- read_csv("data/derived/final_table.csv")

# Keep one row per gene (IMPORTANT)
final_genes <- final_table %>%
  distinct(Ssal, .keep_all = TRUE)

# Remove incomplete rows
final_genes <- final_genes %>%
  filter(!is.na(comp_gill), !is.na(comp_liver))

# Define switching (YOUR logic, unchanged)
final_genes <- final_genes %>%
  mutate(
    switch_status = case_when(
      comp_gill == "A" & comp_liver == "A" ~ "A→A",
      comp_gill == "B" & comp_liver == "B" ~ "B→B",
      comp_gill == "A" & comp_liver == "B" ~ "A→B",
      comp_gill == "B" & comp_liver == "A" ~ "B→A"
    )
  )

# Check counts 
table(final_genes$switch_status)

# Save switching counts 
switch_counts <- final_genes %>%
  count(switch_status)

write_csv(
  switch_counts,
  "data/derived/compartment_switch_counts.csv")

