# 03_compartment_distributions.R
# Plot A/B compartment distributions in gill and liver

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# Load data

final_table <- read_csv("data/derived/final_table.csv")

# Keep one row per gene (CRITICAL)
final_genes <- final_table %>%
  distinct(Ssal, .keep_all = TRUE)

# Panel A: Global A/B distribution

pn <- ggplot(
  final_genes %>%
    pivot_longer(
      cols = c(comp_gill, comp_liver),
      names_to = "tissue",
      values_to = "compartment"
    ),
  aes(x = tissue, fill = compartment)
) +
  scale_fill_manual(values = c("A" = "#BF2C34", "B" = "#799FCB")) +
  geom_bar(position = "dodge") +
  scale_x_discrete(
    labels = c(
      comp_gill = "Gill",
      comp_liver = "Liver"
    )
  ) +
  labs(
    title = "A/B Compartment Distribution",
    x = "Tissue",
    y = "Number of Genes",
    fill = "Compartment"
  ) +
  theme_bw()

# Panel B: Chromosome-level A/B

chrom_long <- final_genes %>%
  pivot_longer(
    cols = c(comp_gill, comp_liver),
    names_to = "tissue",
    values_to = "compartment"
  ) %>%
  filter(!is.na(compartment)) %>%
  mutate(
    chrom = factor(chrom, levels = paste0("chr", 1:29)),
    tissue = factor(
      tissue,
      levels = c("comp_gill", "comp_liver"),
      labels = c("Gill", "Liver")
    )
  )

chrom_counts <- chrom_long %>%
  group_by(chrom, tissue, compartment) %>%
  summarise(n = n(), .groups = "drop")

cols_ab <- c("A" = "#BF2C34", "B" = "#799FCB")

plot_comp <- function(df, title_text) {
  ggplot(df,
         aes(x = chrom, y = n, fill = compartment)) +
    geom_col(
      position = position_dodge(width = 0.8),
      width = 0.8,
      color = "grey20",
      linewidth = 0.2
    ) +
    scale_fill_manual(values = cols_ab) +
    labs(
      title = title_text,
      x = "Chromosome",
      y = "Number of Compartments",
      fill = "Compartment"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

p_gill_genes <- chrom_counts %>%
  filter(tissue == "Gill") %>%
  plot_comp("Chromosome level A/B Distribution — Gill")

p_liver_genes <- chrom_counts %>%
  filter(tissue == "Liver") %>%
  plot_comp("Chromosome level A/B Distribution — Liver")

# Combine panels

comp_dist <- (pn | (p_gill_genes / p_liver_genes)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0, 1)
  )

# Save figure

ggsave(
  filename = "figures/Figure_AB_compartments.tiff",
  plot = res1,
  width = 24,
  height = 14,
  units = "cm",
  dpi = 300,
  compression = "lzw"
)
