# Title: Multiomics for Cytoscape analysis
# Author: Dr. Alice M. Godden

# ============================================
# MASTER MULTI-OMIC NETWORK PIPELINE 
# ============================================

library(dplyr)
library(readr)
library(stringr)
library(AnnotationDbi)
library(org.Dr.eg.db)

# ============================================
# 1. ---------- EDGES ----------
# ============================================

# ---------- miRANDA (miRNA → gene SYMBOL) ----------
miranda_raw <- read.table(
  "output_tab_with_csv_TESTES.txt",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

miranda_edges <- data.frame(
  source = miranda_raw$Seq1,
  target = mapIds(
    org.Dr.eg.db,
    keys = miranda_raw$Parent,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  ),
  score = miranda_raw$Tot.Score
) %>%
  dplyr::filter(score > 150, !is.na(target)) %>%
  dplyr::distinct(source, target) %>%
  dplyr::mutate(interaction = "miRNA_gene")

# ---------- FishPi (piRNA → TE family) ----------


fishpi <- read.csv("complementary_TE_list_HS_te_allpiRNAs.csv")

# ✅ Clean names just in case

colnames(fishpi) 

# force correct column name
colnames(fishpi)[grepl("weight", colnames(fishpi))] <- "weight"
# ✅ Build edges (now simple)

fishpi_edges <- fishpi %>%
  dplyr::select(
    piRNA.Name,
    TE.Name,
    weight
  ) %>%
  dplyr::rename(
    source = piRNA.Name,
    TE_raw = TE.Name
  ) %>%
  dplyr::mutate(
    weight = as.numeric(weight),
    TE_raw = gsub("^>", "", TE_raw),
    TE_family = stringr::str_split_fixed(TE_raw, "::", 2)[,1]
  ) %>%
  dplyr::select(source, target = TE_family, weight) %>%
  dplyr::distinct(source, target, weight) %>%
  dplyr::filter(!is.na(weight), weight >= 20) %>%
  dplyr::mutate(interaction = "piRNA_TE")

# ---------- FishTEA (TE → gene SYMBOL) ----------
fishtea <- read.csv("fishtea_ov_matched_use.csv")

fishtea_edges <- fishtea %>%
  dplyr::select(TE_name, Gene_name, TE_class) %>%
  dplyr::rename(source = TE_name, target = Gene_name) %>%
  dplyr::filter(target != "", !is.na(target),
                TE_class %in% c("LINE","LTR","DNA")) %>%
  dplyr::mutate(
    source = str_split_fixed(source, ":", 2)[,1]   # ✅ match FishPi TE family
  ) %>%
  dplyr::distinct(source, target) %>%
  dplyr::mutate(interaction = "TE_gene")

# ---------- COMBINE EDGES ----------
all_edges <- bind_rows(
  miranda_edges,
  fishpi_edges,
  #fishtea_edges
)

# ============================================
# 2. ---------- NODES ----------
# ============================================

clean_nodes <- function(df, type, id_col, fc_col) {
  
  colnames(df) <- trimws(colnames(df))
  
  df %>%
    dplyr::rename(feature_id = !!id_col,
                  log2FC = !!fc_col) %>%
    dplyr::mutate(
      feature_type = type,
      significant = padj < 0.05 & abs(log2FC) > 1,
      direction = case_when(
        log2FC > 1 ~ "up",
        log2FC < -1 ~ "down",
        TRUE ~ "neutral"
      )
    ) %>%
    dplyr::filter(significant)
}

# ---------- piRNA ----------
pirna_nodes <- clean_nodes(
  read_csv("12-male_tesmall_pirna_volcano.csv"),
  "piRNA", rlang::sym("piRNA"), rlang::sym("log2FoldChange")
)

# ---------- miRNA ----------
mirna_nodes <- clean_nodes(
  read_tsv("8-TESTES_DEDUP_MATURE_TVsC_NA_use2.tsv"),
  "miRNA", rlang::sym("miRNA"), rlang::sym("log2FoldChange")
)

# ---------- TE ----------
te_raw <- read_csv("4-te_tetrans_deseq2.csv")
colnames(te_raw)[1] <- "feature_id"

te_nodes <- te_raw %>%
  dplyr::rename(log2FC = log2FoldChange) %>%
  dplyr::filter(baseMean > 10) %>%
  dplyr::mutate(
    feature_type = "TE",
    feature_id = str_split(feature_id, ":", simplify = TRUE)[,1],  # ✅ FIX
    significant = padj < 0.05 & abs(log2FC) > 1,
    direction = case_when(
      log2FC > 1 ~ "up",
      log2FC < -1 ~ "down",
      TRUE ~ "neutral"
    )
  ) %>%
  dplyr::filter(significant)

# ---------- genes (ENSEMBL → SYMBOL) ----------
gene_raw <- read_csv("2-te_tocontrol_deseq2.csv")

colnames(gene_raw)[1] <- "feature_id"

symbols <- mapIds(
  org.Dr.eg.db,
  keys = gene_raw$feature_id,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_nodes <- gene_raw %>%
  dplyr::mutate(feature_id = symbols) %>%
  dplyr::filter(!is.na(feature_id)) %>%
  dplyr::rename(log2FC = log2FoldChange) %>%
  dplyr::filter(baseMean > 10) %>%
  dplyr::mutate(
    feature_type = "gene",
    significant = padj < 0.05 & abs(log2FC) > 1,
    direction = case_when(
      log2FC > 1 ~ "up",
      log2FC < -1 ~ "down",
      TRUE ~ "neutral"
    )
  ) %>%
  dplyr::filter(significant)

# ---------- COMBINE NODES ----------
all_nodes <- bind_rows(
  gene_nodes,
  mirna_nodes,
  pirna_nodes,
  te_nodes
)

# ============================================
# 3. ---------- FILTER EDGES ----------
# ============================================

valid_nodes <- all_nodes$feature_id

all_edges <- all_edges %>%
  dplyr::filter(
    source %in% valid_nodes &
      target %in% valid_nodes
  )


# ✅ DO NOT OVER-FILTER HERE

# ============================================
# 4. ---------- SAVE ----------
# ============================================

write_csv(all_nodes, "cleaned/male_FINAL_nodes.csv")
write_csv(all_edges, "cleaned/male_FINAL_edges.csv")

cat("Nodes:", nrow(all_nodes), "\n")
cat("Edges:", nrow(all_edges), "\n")

# find key regulators
library(dplyr)

hub_nodes <- all_edges %>%
  group_by(source) %>%
  summarise(degree = n()) %>%
  arrange(desc(degree))

head(hub_nodes, 20)

# graph
library(igraph)

g <- graph_from_data_frame(all_edges, vertices = all_nodes)

plot(g, vertex.size = 5, vertex.label = NA)

