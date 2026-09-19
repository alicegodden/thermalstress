library(EnhancedVolcano)
library(viridis)
library(dplyr)
library(ggrepel)

deseq_results <- read_tsv('VOLCANO_annotated_te_tocontrol_deseq2.tsv')
highlight_genes <- c("hspa8b", "hsp90aa1.2", "hsp70.3", "marcksl1a")

# top 30 by padj, but EXCLUDE the highlight genes so they aren't double-labelled
genes_to_label <- deseq_results %>% arrange(padj) %>%
  pull(external_gene_name) %>% head(50)
genes_to_label <- setdiff(genes_to_label, highlight_genes)   # <-- omit highlights here

p <- EnhancedVolcano(
  deseq_results,
  lab = deseq_results$external_gene_name,
  x = 'log2FoldChange', y = 'padj',
  selectLab = genes_to_label,
  col = rocket(6),
  pCutoff = 50e-03, FCcutoff = log2(1.5),
  ylim = c(0, 10), xlim = c(-10, 30),
  min.segment.length = 0.1,
  labFace = "bold", labSize = 2.5, pointSize = 1.5,
  drawConnectors = TRUE, widthConnectors = 0.2,
  title = 'Enhanced Volcano',
  subtitle = 'RNA-seq: Testes thermal stress'
)

# overlay ONLY the four highlight genes, in larger bold font
hl <- deseq_results %>% filter(external_gene_name %in% highlight_genes)
p + geom_text_repel(
  data = hl,
  aes(x = log2FoldChange, y = -log10(padj), label = external_gene_name),
  size = 3.25, fontface = "bold", colour = "royalblue4",
  min.segment.length = 0.1, box.padding = 0.6, max.overlaps = Inf
)
