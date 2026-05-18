# Title: Regulon and transcription factor analysis
# Author: Dr. Alice M. Godden

# =========================
# Install + load packages
# =========================
install.packages("remotes")

library(pheatmap)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(dorothea)
library(viper)
# =========================
# 1. Load expression data
# =========================
expr <- read.csv("normalised_counts_genes_male.csv", row.names = 1)
expr <- as.matrix(expr)

# =========================
# 2. Convert ENS → SYMBOL
# =========================
symbols <- mapIds(
  org.Dr.eg.db,
  keys = rownames(expr),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

valid <- !is.na(symbols)
expr <- expr[valid, ]
symbols <- symbols[valid]

symbols <- make.unique(symbols)
rownames(expr) <- symbols

# =========================
# 3. Load ZEBRAFISH TF list ✅ (IMPROVED)
# =========================

zf_tfs <- read.table(
  "Danio_rerio_TF.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Extract Ensembl IDs and Symbols
tf_ensembl <- zf_tfs$Ensembl
tf_symbols <- zf_tfs$Symbol

# =========================
# Match using ENSEMBL first (best)
# =========================

# Keep genes where ENSEMBL matches
tf_expr <- expr[rownames(expr) %in% tf_ensembl, ]

cat("TFs matched by ENSEMBL:", nrow(tf_expr), "\n")

# =========================
# OPTIONAL: fallback using SYMBOL (adds extras)
# =========================

extra_tfs <- intersect(tf_symbols, rownames(expr))

tf_expr_extra <- expr[rownames(expr) %in% extra_tfs, ]

# Combine (avoid duplicates)
tf_expr <- rbind(tf_expr, tf_expr_extra[!rownames(tf_expr_extra) %in% rownames(tf_expr), ])

cat("Total TFs after combining:", nrow(tf_expr), "\n")

# =========================
# 5. Clean matrix
# =========================
tf_expr <- tf_expr[
  apply(tf_expr, 1, function(x) all(is.finite(x)) && var(x) > 0),
]

# =========================
# 6. Define groups
# =========================
control_idx <- grep("^MC", colnames(expr))
temp_idx <- grep("^MT", colnames(expr))

annotation_col <- data.frame(
  Temperature = ifelse(grepl("^MC", colnames(tf_expr)), "Control", "Temp")
)

rownames(annotation_col) <- colnames(tf_expr)
annotation_col$Temperature <- factor(annotation_col$Temperature)

# =========================
# 7. Order samples
# =========================
order_idx <- order(annotation_col$Temperature)

tf_expr <- tf_expr[, order_idx]
annotation_col <- annotation_col[order_idx, , drop = FALSE]

# =========================
# 8. Compute logFC
# =========================
mean_control <- rowMeans(expr[, control_idx])
mean_temp <- rowMeans(expr[, temp_idx])

logFC <- log2(mean_temp + 1) - log2(mean_control + 1)
names(logFC) <- rownames(expr)

# =========================
# 9. Statistical testing
# =========================
pvals <- apply(tf_expr, 1, function(gene) {
  t.test(gene[temp_idx], gene[control_idx])$p.value
})

padj <- p.adjust(pvals, method = "BH")

tf_stats <- data.frame(
  Gene = rownames(tf_expr),
  logFC = logFC[rownames(tf_expr)],
  p_value = pvals,
  adj_p_value = padj
)

# Label significance
tf_stats$significance <- "Not Sig"
tf_stats$significance[tf_stats$adj_p_value < 0.05 & tf_stats$logFC > 0] <- "Up"
tf_stats$significance[tf_stats$adj_p_value < 0.05 & tf_stats$logFC < 0] <- "Down"

# =========================
# 10. Save CSV ✅
# =========================
output_file <- "Zebrafish_TF_DE_results.csv"
write.csv(tf_stats, output_file, row.names = FALSE)

cat("CSV saved to:", normalizePath(output_file), "\n")

# =========================
# 11. Select TOP 20 TFs ✅
# =========================
tf_stats_sorted <- tf_stats[order(tf_stats$adj_p_value, -abs(tf_stats$logFC)), ]

sig_tfs <- tf_stats_sorted[tf_stats_sorted$adj_p_value < 0.05, ]

top_tf_genes <- sig_tfs$Gene

# Fill to 20 if needed
if (length(top_tf_genes) < 20) {
  remaining <- tf_stats_sorted$Gene[!tf_stats_sorted$Gene %in% top_tf_genes]
  needed <- 20 - length(top_tf_genes)
  top_tf_genes <- c(top_tf_genes, head(remaining, needed))
}

top_tf_genes <- head(top_tf_genes, 20)

cat("TFs in heatmap:", length(top_tf_genes), "\n")

# =========================
# 12. Plot ALL TFs
# =========================
colors <- colorRampPalette(c("navy","white","firebrick3"))(100)

pheatmap(
  tf_expr,
  scale = "row",
  annotation_col = annotation_col,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = FALSE,
  gaps_col = sum(annotation_col$Temperature == "Control"),
  color = colors,
  main = "All zebrafish transcription factors"
)

# =========================
# 13. Plot TOP 20 TFs ✅
# =========================
tf_expr_top <- tf_expr[top_tf_genes, , drop = FALSE]

# Order by logFC
tf_expr_top <- tf_expr_top[order(logFC[rownames(tf_expr_top)], decreasing = TRUE), ]

annotation_col_top <- annotation_col[colnames(tf_expr_top), , drop = FALSE]

# Add significance stars
sig_lookup <- tf_stats$adj_p_value
names(sig_lookup) <- tf_stats$Gene

stars <- ifelse(sig_lookup[rownames(tf_expr_top)] < 0.05, "*", "")
rownames(tf_expr_top) <- paste0(rownames(tf_expr_top), stars)

pheatmap(
  tf_expr_top,
  scale = "row",
  annotation_col = annotation_col_top,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  gaps_col = sum(annotation_col_top$Temperature == "Control"),
  color = colors,
  show_rownames = TRUE,
  main = "Top 20 zebrafish TFs (* = FDR < 0.05)"
)

# =========================
# 14. Summary
# =========================
cat("Significant TFs (FDR < 0.05):", sum(tf_stats$adj_p_value < 0.05), "\n")


        # no significant results were found
