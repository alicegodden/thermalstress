# =====================================================================
# Enrichment of heat-induced (temperature-unique) de novo TE insertions
# in open chromatin, ovary vs testis.
# Insertions: RetroSeq treatment-unique VCFs (control already excluded)
# Open chromatin: DANIO-CODE DOPES open-region BED
# =====================================================================
library(dplyr); library(tidyr); library(ggplot2)
library(GenomicRanges)

GENOME_BP <- 1.37e9   # zebrafish GRCz11 total assembly length (~1.37 Gb)

# ---- read RetroSeq insertion VCF -> chrom, pos, TE family -----------
read_ins <- function(path, tissue) {
  v <- readLines(path)
  v <- v[!grepl("^#", v)]
  parts <- strsplit(v, "\t")
  chrom <- sapply(parts, `[`, 1)
  pos   <- as.integer(sapply(parts, `[`, 2))
  info  <- sapply(parts, `[`, 8)
  fam   <- sub(".*MEINFO=([^,]+),.*", "\\1", info)
  data.frame(chrom = sub("^chr", "", chrom), pos = pos,
             family = fam, Tissue = tissue, stringsAsFactors = FALSE)
}

# ---- read open-region BED ------------------------------------------
read_open <- function(path) {
  b <- read.table(path, header = FALSE, sep = "\t", comment.char = "",
                  stringsAsFactors = FALSE)
  # header line begins with #chrom; drop if present
  if (grepl("chrom", b[1,1], ignore.case = TRUE)) b <- b[-1, ]
  data.frame(chrom = sub("^chr", "", b[[1]]),
             start = as.integer(b[[2]]), end = as.integer(b[[3]]))
}

open_bed <- read_open("daniocode_hub_280355_dopes_all.txt")
open_gr  <- GRanges(open_bed$chrom, IRanges(open_bed$start, open_bed$end))
open_frac <- sum(width(reduce(open_gr))) / GENOME_BP   # fraction of genome that is open

ins <- bind_rows(
  read_ins("male_exp_unique_8_win100_filterpy_gq750_fl8.vcf",   "Testis"),
  read_ins("female_exp_unique_8_win100_filterpy_gq750_fl8.vcf", "Ovary"))

# ---- classify each insertion open / closed by overlap --------------
ins_gr <- GRanges(ins$chrom, IRanges(ins$pos, ins$pos))
ins$Chromatin <- ifelse(overlapsAny(ins_gr, open_gr), "Open", "Closed")

# =====================================================================
# ENRICHMENT TEST — observed open fraction vs genome coverage, per sex
# Binomial test: are insertions in open chromatin MORE than expected
# from the genome-wide coverage of open regions?
# =====================================================================
enrich <- ins %>% group_by(Tissue) %>%
  summarise(N = n(), N_open = sum(Chromatin == "Open"), .groups = "drop") %>%
  rowwise() %>%
  mutate(pct_open = 100 * N_open / N,
         expected_pct = 100 * open_frac,
         enrichment = (N_open / N) / open_frac,
         p_value = binom.test(N_open, N, open_frac, alternative = "greater")$p.value,
         ci_lo = 100 * binom.test(N_open, N, open_frac)$conf.int[1],
         ci_hi = 100 * binom.test(N_open, N, open_frac)$conf.int[2]) %>%
  ungroup()

cat(sprintf("\nOpen chromatin covers %.3f%% of the genome (DANIO-CODE DOPES).\n\n",
            100 * open_frac))
print(as.data.frame(enrich))

# =====================================================================
# FIGURE: observed % open per sex vs expected (dashed line), fold + p
# =====================================================================
p <- ggplot(enrich, aes(Tissue, pct_open, fill = Tissue)) +
  geom_col(width = .6, alpha = .9) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .15, linewidth = .8) +
  geom_hline(yintercept = 100 * open_frac, linetype = "dashed", color = "grey30") +
  annotate("text", x = 0.5, y = 100 * open_frac + 0.12, hjust = 0,
           label = sprintf("expected by coverage (%.2f%%)", 100 * open_frac),
           size = 3.4, fontface = "bold", color = "grey30") +
  geom_text(aes(y = ci_hi + 0.15,
                label = sprintf("%.1f\u00d7 enriched\nP = %.2g", enrichment, p_value)),
            fontface = "bold", size = 4, lineheight = .9) +
  scale_fill_manual(values = c("Ovary" = "#d98f52", "Testis" = "#c85a28"), guide = "none") +
  labs(x = NULL, y = "Insertions in open chromatin (%)",
       title = "Heat-induced TE insertions are enriched in open chromatin",
       subtitle = "Temperature-unique de novo insertions vs genome-wide open-region coverage") +
  theme_classic() +
  theme(plot.title = element_text(size = 15, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text  = element_text(size = 12, face = "bold"))
ggsave("fig_TE_open_enrichment.png", p, width = 7, height = 6, dpi = 600)

write.csv(enrich, "TE_open_chromatin_enrichment.csv", row.names = FALSE)
message("Done: fig_TE_open_enrichment.png, TE_open_chromatin_enrichment.csv")



# =====================================================================
# Enrichment of heat-induced (temperature-unique) de novo TE insertions
# in open chromatin, ovary vs testis.
# Insertions: RetroSeq treatment-unique VCFs (control already excluded)
# Open chromatin: DANIO-CODE DOPES open-region BED
# =====================================================================
library(dplyr); library(tidyr); library(ggplot2)
library(GenomicRanges)

GENOME_BP <- 1.37e9   # zebrafish GRCz11 total assembly length (~1.37 Gb)

# ---- read RetroSeq insertion VCF -> chrom, pos, TE family -----------
read_ins <- function(path, tissue) {
  v <- readLines(path)
  v <- v[!grepl("^#", v)]
  parts <- strsplit(v, "\t")
  chrom <- sapply(parts, `[`, 1)
  pos   <- as.integer(sapply(parts, `[`, 2))
  info  <- sapply(parts, `[`, 8)
  fam   <- sub(".*MEINFO=([^,]+),.*", "\\1", info)
  data.frame(chrom = sub("^chr", "", chrom), pos = pos,
             family = fam, Tissue = tissue, stringsAsFactors = FALSE)
}

# ---- read open-region BED ------------------------------------------
read_open <- function(path) {
  b <- read.table(path, header = FALSE, sep = "\t", comment.char = "",
                  stringsAsFactors = FALSE)
  # header line begins with #chrom; drop if present
  if (grepl("chrom", b[1,1], ignore.case = TRUE)) b <- b[-1, ]
  data.frame(chrom = sub("^chr", "", b[[1]]),
             start = as.integer(b[[2]]), end = as.integer(b[[3]]))
}

open_bed <- read_open("daniocode_hub_280355_dopes_all.txt")
open_gr  <- GRanges(open_bed$chrom, IRanges(open_bed$start, open_bed$end))
open_frac <- sum(width(reduce(open_gr))) / GENOME_BP   # fraction of genome that is open

ins <- bind_rows(
  read_ins("male_exp_unique_8_win100_filterpy_gq750_fl8.vcf",   "Testis"),
  read_ins("1789584653432_female_exp_unique_8_win100_filterpy_gq750_fl8.vcf", "Ovary"))

# ---- classify each insertion open / closed by overlap --------------
ins_gr <- GRanges(ins$chrom, IRanges(ins$pos, ins$pos))
ins$Chromatin <- ifelse(overlapsAny(ins_gr, open_gr), "Open", "Closed")

# =====================================================================
# ENRICHMENT TEST — observed open fraction vs genome coverage, per sex
# Binomial test: are insertions in open chromatin MORE than expected
# from the genome-wide coverage of open regions?
# =====================================================================
enrich <- ins %>% group_by(Tissue) %>%
  summarise(N = n(), N_open = sum(Chromatin == "Open"), .groups = "drop") %>%
  rowwise() %>%
  mutate(pct_open = 100 * N_open / N,
         expected_pct = 100 * open_frac,
         enrichment = (N_open / N) / open_frac,
         p_value = binom.test(N_open, N, open_frac, alternative = "greater")$p.value,
         ci_lo = 100 * binom.test(N_open, N, open_frac)$conf.int[1],
         ci_hi = 100 * binom.test(N_open, N, open_frac)$conf.int[2]) %>%
  ungroup()

cat(sprintf("\nOpen chromatin covers %.3f%% of the genome (DANIO-CODE DOPES).\n\n",
            100 * open_frac))
print(as.data.frame(enrich))


# =====================================================================
# COUNTS TABLE — open vs closed insertions, both sexes (easy to read)
# =====================================================================
counts_tab <- ins %>%
  group_by(Tissue, Chromatin) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Chromatin, values_from = n, values_fill = 0) %>%
  mutate(Total       = Open + Closed,
         pct_Open    = round(100 * Open   / Total, 2),
         pct_Closed  = round(100 * Closed / Total, 2),
         expected_pct_Open = round(100 * open_frac, 2),
         fold_Open   = round((Open / Total) / open_frac, 2)) %>%
  select(Tissue, Total, Open, Closed, pct_Open, pct_Closed,
         expected_pct_Open, fold_Open)

cat("\n===== Open vs closed insertion counts =====\n")
print(as.data.frame(counts_tab))
write.csv(counts_tab, "TE_open_closed_counts.csv", row.names = FALSE)

# =====================================================================
# FIGURE: observed % open per sex vs expected (dashed line), fold + p
# =====================================================================
p <- ggplot(enrich, aes(Tissue, pct_open, fill = Tissue)) +
  geom_col(width = .6, alpha = .9) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .15, linewidth = .8) +
  geom_hline(yintercept = 100 * open_frac, linetype = "dashed", color = "grey30") +
  annotate("text", x = 0.5, y = 100 * open_frac + 0.12, hjust = 0,
           label = sprintf("expected by coverage (%.2f%%)", 100 * open_frac),
           size = 3.4, fontface = "bold", color = "grey30") +
  geom_text(aes(y = ci_hi + 0.15,
                label = sprintf("%.1f\u00d7 enriched\nP = %.2g", enrichment, p_value)),
            fontface = "bold", size = 4, lineheight = .9) +
  scale_fill_manual(values = c("Ovary" = "#d98f52", "Testis" = "#c85a28"), guide = "none") +
  labs(x = NULL, y = "Insertions in open chromatin (%)",
       title = "Heat-induced TE insertions are enriched in open chromatin",
       subtitle = "Temperature-unique de novo insertions vs genome-wide open-region coverage") +
  theme_classic() +
  theme(plot.title = element_text(size = 15, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text  = element_text(size = 12, face = "bold"))
ggsave("fig_TE_open_enrichment.png", p, width = 7, height = 6, dpi = 600)


# =====================================================================
# FIGURE 2: open AND closed counts side by side (stacked bar), per sex
# =====================================================================
long_ct <- counts_tab %>%
  select(Tissue, Open, Closed) %>%
  pivot_longer(c(Open, Closed), names_to = "Chromatin", values_to = "n") %>%
  mutate(Chromatin = factor(Chromatin, levels = c("Closed","Open")))

p2 <- ggplot(long_ct, aes(Tissue, n, fill = Chromatin)) +
  geom_col(width = .6, alpha = .9) +
  geom_text(aes(label = n), position = position_stack(vjust = .5),
            fontface = "bold", size = 4, color = "white") +
  scale_fill_manual(values = c("Closed" = "#6a7488", "Open" = "#e8b64c")) +
  labs(x = NULL, y = "Temperature-unique de novo insertions",
       title = "De novo TE insertions by chromatin state and sex",
       subtitle = "Open (gold) vs closed (grey); open chromatin = 0.39% of the genome") +
  theme_classic() +
  theme(plot.title = element_text(size = 15, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "grey40"),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text  = element_text(size = 12, face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"))
ggsave("fig_TE_open_closed_counts.png", p2, width = 7, height = 6, dpi = 600)

write.csv(enrich, "TE_open_chromatin_enrichment.csv", row.names = FALSE)
message("Done: fig_TE_open_enrichment.png, fig_TE_open_closed_counts.png, TE_open_chromatin_enrichment.csv, TE_open_closed_counts.csv")

