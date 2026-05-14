# 07_ohnolog_analysis.R
# Ohnolog analysis: classification, pairs, regulatory divergence, GO

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
})

# 1. LOAD DATA

final_table <- read_csv("data/derived/final_table.csv")

final_genes <- final_table %>%
  distinct(Ssal, .keep_all = TRUE)

ohnolog_pairs_ens <- read_csv("data/derived/ohnolog_pairs_ens.csv")
ohnolog_ids       <- read_csv("data/derived/ohnolog_ids.csv")
singleton_ids     <- read_csv("data/derived/singleton_ids.csv")

# 2. ASSIGN OHNOLOG STATUS (GENE LEVEL)

final_genes <- final_genes %>%
  left_join(ohnolog_ids, by = c("Ssal" = "ensembl_gene_id")) %>%
  mutate(is_ohnolog = !is.na(ensembl_gene_id)) %>%
  select(-ensembl_gene_id) %>%
  left_join(singleton_ids, by = c("Ssal" = "ensembl_gene_id")) %>%
  mutate(is_singleton = !is.na(ensembl_gene_id)) %>%
  select(-ensembl_gene_id) %>%
  mutate(
    ohnolog_status = case_when(
      is_ohnolog   ~ "ohnolog",
      is_singleton ~ "singleton",
      dup_type == "duplicate" ~ "other_duplicate",
      TRUE ~ "unknown"
    )
  )

# 3. BUILD CLEAN 1:1 OHNOLOG PAIRS

pairs_raw <- ohnolog_pairs_ens %>%
  transmute(
    gene1 = str_squish(ensembl_A),
    gene2 = str_squish(ensembl_B)
  ) %>%
  filter(!is.na(gene1), !is.na(gene2), gene1 != gene2)

pairs_clean <- pairs_raw %>%
  mutate(
    gmin = pmin(gene1, gene2),
    gmax = pmax(gene1, gene2)
  ) %>%
  distinct(gmin, gmax) %>%
  rename(gene1 = gmin, gene2 = gmax)

pairs_1to1 <- pairs_clean %>%
  group_by(gene1) %>% filter(n() == 1) %>% ungroup() %>%
  group_by(gene2) %>% filter(n() == 1) %>% ungroup()

# 4. ADD COMPARTMENTS TO PAIRS

pairs_comp <- pairs_1to1 %>%
  left_join(final_genes %>% select(Ssal, comp_gill, comp_liver),
            by = c("gene1" = "Ssal")) %>%
  rename(comp_gill_1 = comp_gill,
         comp_liver_1 = comp_liver) %>%
  left_join(final_genes %>% select(Ssal, comp_gill, comp_liver),
            by = c("gene2" = "Ssal")) %>%
  rename(comp_gill_2 = comp_gill,
         comp_liver_2 = comp_liver)

# 5. CLASSIFY SAME / DIFFERENT COMPARTMENTS

pairs_comp <- pairs_comp %>%
  filter(
    comp_gill_1 %in% c("A","B"),
    comp_gill_2 %in% c("A","B"),
    comp_liver_1 %in% c("A","B"),
    comp_liver_2 %in% c("A","B")
  ) %>%
  mutate(
    gill_status  = ifelse(comp_gill_1 == comp_gill_2, "Same", "Different"),
    liver_status = ifelse(comp_liver_1 == comp_liver_2, "Same", "Different")
  )

# 6. ADD THETA DATA

theta_map <- final_genes %>%
  filter(dup_type == "duplicate") %>%
  select(Ssal, thetaShift, isSig)

pairs_theta <- pairs_1to1 %>%
  left_join(theta_map, by = c("gene1" = "Ssal")) %>%
  rename(theta_1 = thetaShift, sig_1 = isSig) %>%
  left_join(theta_map, by = c("gene2" = "Ssal")) %>%
  rename(theta_2 = thetaShift, sig_2 = isSig)

pairs_theta <- pairs_theta %>%
  mutate(
    theta_pair = pmax(abs(theta_1), abs(theta_2), na.rm = TRUE),
    sig_pair   = sig_1 | sig_2
  )

# 7. COMBINE PAIRS + COMPARTMENTS + THETA

pairs_final <- pairs_comp %>%
  left_join(pairs_theta, by = c("gene1", "gene2"))

# 8. BUILD FINAL ANALYSIS TABLE

theta_df <- bind_rows(
  data.frame(
    tissue = "Gill",
    comp_class = pairs_final$gill_status,
    theta = pairs_final$theta_pair,
    sig = pairs_final$sig_pair
  ),
  data.frame(
    tissue = "Liver",
    comp_class = pairs_final$liver_status,
    theta = pairs_final$theta_pair,
    sig = pairs_final$sig_pair
  )
) %>%
  filter(!is.na(comp_class))

# 9. STATISTICS

# Fisher test
fisher_gill  <- fisher.test(table(pairs_comp$gill_status))
fisher_liver <- fisher.test(table(pairs_comp$liver_status))

# Wilcoxon (explicit)
gill_same <- theta_df$theta[theta_df$tissue=="Gill" & theta_df$comp_class=="Same"]
gill_diff <- theta_df$theta[theta_df$tissue=="Gill" & theta_df$comp_class=="Different"]

liver_same <- theta_df$theta[theta_df$tissue=="Liver" & theta_df$comp_class=="Same"]
liver_diff <- theta_df$theta[theta_df$tissue=="Liver" & theta_df$comp_class=="Different"]

wilcox_gill  <- wilcox.test(gill_same, gill_diff)
wilcox_liver <- wilcox.test(liver_same, liver_diff)

# 10. GO ENRICHMENT

ann <- readRDS("data/derived/ssal_GO_BP_biomart.rds")

ann <- ann %>%
  filter(go_id != "", !is.na(go_id))

TERM2GENE <- ann %>%
  select(go_id, ensembl_gene_id)

TERM2NAME <- ann %>%
  select(go_id, name_1006) %>%
  distinct()

ohnolog_genes <- final_genes %>%
  filter(ohnolog_status == "ohnolog") %>%
  pull(Ssal) %>%
  unique()

ego <- enricher(
  gene = ohnolog_genes,
  universe = unique(TERM2GENE$ensembl_gene_id),
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME
)