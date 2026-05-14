# 06_ohnolog_mapping.R
# Build LOC → Ensembl mapping and define ohnolog / singleton sets

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# Load annotations

ncbi_gff <- import("data/raw/GCF_905237065.1_Ssal_v3.1_genomic.gff")
ens_gff  <- import("data/raw/Salmo_salar.Ssal_v3.1.106.gff3")

# NCBI genes (LOC)
ncbi_genes <- ncbi_gff[ncbi_gff$type == "gene"]

gr_loc <- GRanges(
  seqnames = seqnames(ncbi_genes),
  ranges   = ranges(ncbi_genes),
  loc_id   = ncbi_genes$Name
)

# Ensembl genes
ens_genes <- ens_gff[ens_gff$type == "gene"]

gr_ens <- GRanges(
  seqnames = seqnames(ens_genes),
  ranges   = ranges(ens_genes),
  ensembl  = ens_genes$gene_id
)

# Rename chromosomes

ncbi_chrs <- paste0("NC_0594", sprintf("%02d", 42:70), ".1")
ens_chrs  <- as.character(1:29)

ncbi_to_ens <- setNames(ens_chrs, ncbi_chrs)

seqlevels(gr_loc) <- plyr::mapvalues(
  seqlevels(gr_loc),
  from = names(ncbi_to_ens),
  to   = ncbi_to_ens,
  warn_missing = FALSE
)

# Overlap + best mapping

hits <- findOverlaps(gr_loc, gr_ens)

ovl_width <- width(
  pintersect(
    gr_loc[queryHits(hits)],
    gr_ens[subjectHits(hits)]
  )
)

map_best <- tibble(
  loc_id  = mcols(gr_loc)$loc_id[queryHits(hits)],
  ensembl = mcols(gr_ens)$ensembl[subjectHits(hits)],
  width   = ovl_width
) %>%
  group_by(loc_id) %>%
  slice_max(width, n = 1, with_ties = FALSE) %>%
  ungroup()

write_csv(map_best, "data/derived/loc_to_ensembl_best.csv")

# Load ohnolog raw data

ohnos_raw <- read.csv("data/raw/Supplementary_Data_13.csv",
                      sep = ";",
                      header = TRUE,
                      skip = 1,
                      fill = TRUE)

# Build ohnolog pairs

pairs_raw <- data.frame(
  locA = ohnos_raw[[3]],
  locB = ohnos_raw[[8]],
  stringsAsFactors = FALSE
)

pairs_raw <- subset(pairs_raw,
                    grepl("^LOC", locA) & grepl("^LOC", locB))

pairs_ens <- merge(pairs_raw, map_best,
                   by.x = "locA", by.y = "loc_id")

pairs_ens <- merge(pairs_ens, map_best,
                   by.x = "locB", by.y = "loc_id",
                   suffixes = c("_A", "_B"))

ohnolog_pairs_ens <- unique(
  pairs_ens[, c("ensembl_A", "ensembl_B")]
)

write_csv(ohnolog_pairs_ens, "data/derived/ohnolog_pairs_ens.csv")

# Ohnolog + singleton IDs

ohnolog_ids <- ohnolog_pairs_ens %>%
  pivot_longer(everything(), values_to = "ensembl_gene_id") %>%
  distinct()

write_csv(ohnolog_ids, "data/derived/ohnolog_ids.csv")

singleton_ids <- ohnos_raw %>%
  filter(!is.na(gene) & is.na(gene.1)) %>%
  transmute(loc_id = gene) %>%
  left_join(map_best, by = "loc_id") %>%
  transmute(ensembl_gene_id = ensembl) %>%
  distinct()

write_csv(singleton_ids, "data/derived/singleton_ids.csv")
