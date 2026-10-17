# =====================================================================
# Age of transcriptionally active vs mobilized TEs under thermal stress
# (reviewer L203-207: plot + stats for transcriptional TE age, and
#  compare with the mobilized subset from RetroSeq)
#
#  - Transcribed set : TElocal/Telescope DESeq2 DE TEs (padj < 0.05), by class
#  - Mobilized set   : RetroSeq private (single-fish) insertions at 34°C,
#                      summed over the 6 heat-exposed fish per organ, by class
#  - Age axis        : RepeatMasker Kimura landscape (divsum), reweighted by
#                      each set's TE-class composition
# Note: curves reflect class composition (the genome-wide age landscape of
# each class, weighted by that class's share of the set), not the ages of the
# individual inserted/transcribed copies.
# Author: Dr Alice M. Godden
# =====================================================================
library(dplyr); library(tidyr); library(ggplot2)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

CLASSES <- c("DNA", "LINE", "LTR", "SINE", "RC", "Satellite")
PADJ    <- 0.05
MOB_LAB <- "Mobilized (private, 34°C)"
TR_LAB  <- "Transcribed (DE)"
BG_LAB  <- "Genome-wide"

# ---- inputs ----------------------------------------------------------
inputs <- list(
  Ovary  = list(de = "11-female_telescope_deseq2.csv", divsum = "FT_all.divsum.csv"),
  Testis = list(de = "10-male_telescope_deseq2.csv",   divsum = "MT_all.divsum.csv"))

# private (single-fish) insertions at 34°C by class, summed over 6 fish
# (from private_TE_by_treatment.tsv)
mob_counts <- list(
  Ovary  = c(DNA = 893,  LINE = 96,  LTR = 193, SINE = 53, RC = 21, Satellite = 15),
  Testis = c(DNA = 1118, LINE = 132, LTR = 244, SINE = 89, RC = 16, Satellite = 39))

# ---- TE family -> broad class (Telescope names are "Family_dupN") ----
map_class <- function(fam) {
  f <- toupper(fam)
  if (grepl("HELITRON", f)) return("RC")
  if (grepl("GYPSY|COPIA|ERV|BEL-|BEL_|PAO|DIRS|LRS|NGARO|GYPSYDR|-LTR", f)) return("LTR")
  if (grepl("CR1|REX|L1-|L2-|L2_|RTE|JOCKEY|I-3|LOOPER|TDR7|TDR23|TDR16|KIRI", f)) return("LINE")
  if (grepl("MOSAT|BRSAT|SAT-|CENSAT|MSAT|\\(", f)) return("Satellite")
  "DNA"   # everything else recognisable is a DNA transposon in zebrafish
}

# ---- transcribed set: DE TE class fractions (padj < PADJ) ------------
de_class_weights <- function(path) {
  d <- read.csv(path, row.names = 1)
  colnames(d) <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  family  <- sub("_dup[0-9]+$", "", rownames(d))
  d$Class <- vapply(family, map_class, character(1))
  sig <- d[!is.na(d$padj) & d$padj < PADJ, ]
  if (nrow(sig) == 0) return(list(w = NULL, n = 0))
  tab <- table(sig$Class)
  list(w = setNames(as.numeric(tab) / sum(tab), names(tab)), n = nrow(sig))
}

mob_weights <- function(organ) { v <- mob_counts[[organ]]; v / sum(v) }

# ---- divsum -> class-weighted age distribution -----------------------
col_to_class <- function(cols) sapply(cols, function(c) {
  if (c == "Div") return(NA)
  for (cl in CLASSES) if (toupper(c) == toupper(cl) ||
                          startsWith(toupper(c), paste0(toupper(cl), "/"))) return(cl)
  NA
})

weighted_age_dist <- function(divsum_path, weights) {
  d <- read.csv(divsum_path, sep = " ", check.names = FALSE)
  d <- d[, !grepl("ARTEFACT|^Unnamed", colnames(d)), drop = FALSE]
  cls <- col_to_class(colnames(d))
  tot <- numeric(nrow(d))
  for (cl in CLASSES) {
    j <- which(cls == cl)
    if (length(j)) tot <- tot + rowSums(as.matrix(d[, j, drop = FALSE])) * (weights[cl] %||% 0)
  }
  tot[tot < 0] <- 0
  data.frame(Div = d$Div, bp = tot)
}

# KS test on age distributions: bp expanded into ~1000 weighted pseudo-samples
# (descriptive only - the p-value depends on this arbitrary sample size)
exp_ages <- function(x) rep(x$Div, times = round(x$bp / max(x$bp) * 1000))

# ---- build per organ + stats -----------------------------------------
plot_df <- list(); stats <- list()
for (organ in names(inputs)) {
  divsum <- inputs[[organ]]$divsum
  mob <- weighted_age_dist(divsum, mob_weights(organ)); mob$set <- MOB_LAB
  bg  <- weighted_age_dist(divsum, setNames(rep(1, length(CLASSES)), CLASSES)); bg$set <- BG_LAB
  de  <- de_class_weights(inputs[[organ]]$de)
  dfs <- list(mob, bg)
  ma_mob <- weighted.mean(mob$Div, mob$bp)
  
  if (!is.null(de$w)) {
    tr <- weighted_age_dist(divsum, de$w); tr$set <- TR_LAB
    dfs <- c(dfs, list(tr))
    ks  <- suppressWarnings(ks.test(exp_ages(tr), exp_ages(mob)))
    ma_tr <- weighted.mean(tr$Div, tr$bp)
    stats[[organ]] <- data.frame(
      Organ = organ, n_DE = de$n, n_mobilized = sum(mob_counts[[organ]]),
      mean_age_transcribed = round(ma_tr, 1), mean_age_mobilized = round(ma_mob, 1),
      KS_D = round(unname(ks$statistic), 3), KS_P = signif(ks$p.value, 3))
  } else {
    stats[[organ]] <- data.frame(
      Organ = organ, n_DE = 0, n_mobilized = sum(mob_counts[[organ]]),
      mean_age_transcribed = NA, mean_age_mobilized = round(ma_mob, 1), KS_D = NA, KS_P = NA)
  }
  d <- bind_rows(dfs); d$Organ <- organ
  d <- d %>% group_by(set) %>% mutate(density = bp / sum(bp)) %>% ungroup()
  plot_df[[organ]] <- d
}

alld     <- bind_rows(plot_df)
stats_df <- bind_rows(stats)
print(stats_df)
write.csv(stats_df, "TE_transcribed_vs_mobilized_age_stats.csv", row.names = FALSE)

alld$set <- factor(alld$set, levels = c(BG_LAB, MOB_LAB, TR_LAB))

# ---- per-facet annotation built from the stats -----------------------
fmt_p <- function(p) ifelse(is.na(p), "NA", ifelse(p < 2.2e-16, "< 2.2e-16", paste("=", signif(p, 2))))
annot <- stats_df %>%
  mutate(label = ifelse(n_DE > 0,
                        sprintf("n = %d DE TEs\nmean age: transcribed %.1f\nvs mobilized %.1f\nKS P %s",
                                n_DE, mean_age_transcribed, mean_age_mobilized, fmt_p(KS_P)),
                        sprintf("no DE TEs (padj < %.2f)\nmean age mobilized %.1f", PADJ, mean_age_mobilized)))
ymax <- max(alld$density[alld$Div <= 40], na.rm = TRUE)

# ---- plot -------------------------------------------------------------
cols <- setNames(c("grey50", "#c85a28", "#2c5f9e"), c(BG_LAB, MOB_LAB, TR_LAB))
p <- ggplot(alld, aes(Div, density, colour = set)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Organ) +
  geom_text(data = annot, aes(x = 24, y = ymax * 0.98, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1,
            fontface = "bold", size = 3, lineheight = 0.95) +
  scale_colour_manual(values = cols, drop = FALSE) +
  coord_cartesian(xlim = c(0, 40)) +
  labs(x = "Kimura substitution level (age)", y = "Relative density", colour = NULL,
       title = "Age of transcriptionally active vs mobilized TEs under thermal stress") +
  theme_classic() +
  theme(plot.title  = element_text(size = 13, face = "bold"),
        axis.title  = element_text(size = 12, face = "bold"),
        axis.text   = element_text(size = 11, face = "bold"),
        strip.text  = element_text(size = 12, face = "bold"),
        legend.text = element_text(face = "bold"),
        legend.position = "bottom")

ggsave("fig_TE_transcribed_vs_mobilized_age.png", p, width = 12, height = 5.5, dpi = 600)
ggsave("fig_TE_transcribed_vs_mobilized_age.pdf", p, width = 12, height = 5.5)
p
