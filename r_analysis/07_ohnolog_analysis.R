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

# =============================
# 2. ASSIGN OHNOLOG STATUS
# =============================

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

write_csv(final_genes, "data/derived/final_genes_ohnolog_status.csv")

# 3. CLEAN OHNOLOG PAIRS (1:1)

pairs_1to1 <- ohnolog_pairs_ens %>%
  transmute(
    gene1 = str_squish(ensembl_A),
    gene2 = str_squish(ensembl_B)
  ) %>%
  filter(!is.na(gene1), !is.na(gene2), gene1 != gene2) %>%
  mutate(
    gmin = pmin(gene1, gene2),
    gmax = pmax(gene1, gene2)
  ) %>%
  distinct(gmin, gmax) %>%
  rename(gene1 = gmin, gene2 = gmax) %>%
  group_by(gene1) %>% filter(n() == 1) %>% ungroup() %>%
  group_by(gene2) %>% filter(n() == 1) %>% ungroup()

write_csv(pairs_1to1, "data/derived/ohnolog_pairs_1to1.csv")

# 4. PAIR COMPARTMENT ANALYSIS

pairs_comp <- pairs_1to1 %>%
  left_join(final_genes %>% select(Ssal, comp_gill, comp_liver),
            by = c("gene1" = "Ssal")) %>%
  rename(comp_gill_1 = comp_gill,
         comp_liver_1 = comp_liver) %>%
  left_join(final_genes %>% select(Ssal, comp_gill, comp_liver),
            by = c("gene2" = "Ssal")) %>%
  rename(comp_gill_2 = comp_gill,
         comp_liver_2 = comp_liver) %>%
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

write_csv(pairs_comp, "data/derived/ohnolog_pairs_compartment.csv")

# 5. REGULATORY DIVERGENCE (THETA)

eve_dup <- final_genes %>%
  filter(dup_type == "duplicate") %>%
  select(Ssal, thetaShift, isSig)

pairs_eve <- pairs_1to1 %>%
  left_join(eve_dup, by = c("gene1" = "Ssal")) %>%
  rename(theta_1 = thetaShift, sig_1 = isSig) %>%
  left_join(eve_dup, by = c("gene2" = "Ssal")) %>%
  rename(theta_2 = thetaShift, sig_2 = isSig) %>%
  mutate(
    theta_pair = pmax(abs(theta_1), abs(theta_2), na.rm = TRUE),
    sig_pair   = sig_1 | sig_2
  )

# Add compartment classification
pairs_eve <- pairs_eve %>%
  left_join(pairs_comp %>% select(gene1, gene2, gill_status, liver_status),
            by = c("gene1", "gene2"))

# Long format
theta_df <- bind_rows(
  tibble(
    tissue = "Gill",
    comp_class = pairs_eve$gill_status,
    theta = pairs_eve$theta_pair,
    sig = pairs_eve$sig_pair
  ),
  tibble(
    tissue = "Liver",
    comp_class = pairs_eve$liver_status,
    theta = pairs_eve$theta_pair,
    sig = pairs_eve$sig_pair
  )
) %>%
  filter(!is.na(comp_class))

write_csv(theta_df, "data/derived/ohnolog_theta_analysis.csv")

# Summary
theta_summary <- theta_df %>%
  group_by(tissue, comp_class) %>%
  summarise(prop_sig = mean(sig), .groups = "drop")

write_csv(theta_summary, "data/derived/ohnolog_theta_summary.csv")

# 6. GO ENRICHMENT

ann <- readRDS("data/derived/ssal_GO_BP_biomart.rds")

ann_clean <- ann %>%
  filter(go_id != "", !is.na(go_id))

TERM2GENE <- ann_clean %>%
  select(go_id, ensembl_gene_id)

TERM2NAME <- ann_clean %>%
  select(go_id, name_1006) %>%
  distinct()

ohnolog_genes <- final_genes %>%
  filter(ohnolog_status == "ohnolog") %>%
  pull(Ssal) %>%
  unique()

background <- unique(TERM2GENE$ensembl_gene_id)

ego <- enricher(
  gene = intersect(ohnolog_genes, background),
  universe = background,
  TERM2GENE = TERM2GENE,
  TERM2NAME = TERM2NAME
)

saveRDS(ego, "data/derived/GO_ohnolog_enrichment.rds")
