# =====================================================================
# Age of transcriptionally active vs mobilized TEs under thermal stress
# (reviewer L203-207: plot + stats for transcriptional TE age, and
#  compare with the mobilized/active subset from RetroSeq)
#
#  - Transcribed set : TETranscripts DESeq2 DE TEs (padj < 0.05), by class
#  - Mobilized set   : RetroSeq temperature-unique insertions, by class
#  - Age axis        : RepeatMasker Kimura landscape (divsum), reweighted
#                      by each activity type's class composition
# Author: Dr Alice M. Godden
# =====================================================================
library(dplyr); library(tidyr); library(ggplot2)

CLASSES <- c("DNA","LINE","LTR","SINE","RC","Satellite")
PADJ <- 0.05

map_class <- function(x) {
  x <- toupper(as.character(x))
  for (cl in CLASSES) if (x == toupper(cl)) return(cl)
  if (grepl("SATELLITE", x)) return("Satellite")
  if (x %in% c("RC","HELITRON")) return("RC")
  NA_character_
}

# ---- transcribed set: DE TE class fractions (padj < 0.05) -----------
de_class_weights <- function(path) {
  d <- read.csv(path, row.names = 1)
  colnames(d) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  cls <- sapply(strsplit(rownames(d), ":"), function(v) map_class(tail(v, 1)))
  d$Class <- cls
  sig <- d[!is.na(d$padj) & d$padj < PADJ & !is.na(d$Class), ]
  if (nrow(sig) == 0) return(list(w = NULL, n = 0))
  tab <- table(sig$Class)
  list(w = setNames(as.numeric(tab)/sum(tab), names(tab)), n = nrow(sig))
}

# ---- mobilized set: RetroSeq class counts (edit to your VCF counts) --
# (broad-class counts of temperature-unique insertions, chr 1-25)
mob_counts <- list(
  Ovary  = c(DNA=524, LINE=25, LTR=85,  SINE=22, RC=5,  Satellite=16),
  Testis = c(DNA=1658,LINE=60, LTR=262, SINE=86, RC=26, Satellite=44))
mob_weights <- function(organ) { v <- mob_counts[[organ]]; v/sum(v) }

# ---- divsum -> weighted age distribution ---------------------------
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
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# ---- build per organ, plot, and KS-test ----------------------------
inputs <- list(
  Ovary  = list(de = "5-_ov_tetrans_deseq2.csv",            divsum = "FT_all.divsum.csv"),
  Testis = list(de = "4-_te_tetrans_deseq2.csv",            divsum = "MT_all.divsum.csv"))

plot_df <- list(); ks_lines <- c()
for (organ in names(inputs)) {
  divsum <- inputs[[organ]]$divsum
  mob <- weighted_age_dist(divsum, mob_weights(organ)); mob$set <- "Mobilized (RetroSeq)"
  bg  <- weighted_age_dist(divsum, setNames(rep(1, length(CLASSES)), CLASSES)); bg$set <- "Genome-wide"
  de  <- de_class_weights(inputs[[organ]]$de)
  dfs <- list(mob, bg)
  if (!is.null(de$w)) {
    tr <- weighted_age_dist(divsum, de$w); tr$set <- sprintf("Transcribed (DE, n=%d)", de$n)
    dfs <- c(dfs, list(tr))
    # KS test on age distributions (expand bp into weighted samples)
    exp_ages <- function(x) rep(x$Div, times = round(x$bp / max(x$bp) * 1000))
    ks <- ks.test(exp_ages(tr), exp_ages(mob))
    ma_tr <- weighted.mean(tr$Div, tr$bp); ma_mob <- weighted.mean(mob$Div, mob$bp)
    ks_lines <- c(ks_lines, sprintf("%s: transcribed mean age %.1f vs mobilized %.1f, KS P = %.3g",
                                    organ, ma_tr, ma_mob, ks$p.value))
  } else {
    ks_lines <- c(ks_lines, sprintf("%s: no DE TEs at padj<%.2f (transcriptional age not computed)", organ, PADJ))
  }
  d <- bind_rows(dfs); d$Organ <- organ
  d <- d %>% group_by(set) %>% mutate(density = bp / sum(bp)) %>% ungroup()
  plot_df[[organ]] <- d
}
alld <- bind_rows(plot_df)
cat("\n===== Transcribed vs mobilized TE age =====\n"); cat(paste(ks_lines, collapse = "\n"), "\n")

p <- ggplot(alld, aes(Div, density, colour = set)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Organ) +
  scale_colour_manual(values = c("Mobilized (RetroSeq)" = "#c85a28",
                                 "Genome-wide" = "grey50"),
                      # transcribed level name varies with n; add by pattern
                      limits = function(x) x) +
  coord_cartesian(xlim = c(0, 40)) +
  labs(x = "Kimura substitution level (age)", y = "Relative density", colour = NULL,
       title = "Age of transcriptionally active vs mobilized TEs under thermal stress") +
  theme_classic() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text  = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(face = "bold"))
ggsave("fig_TE_transcribed_vs_mobilized_age.png", p, width = 12, height = 5.5, dpi = 600)
message("Done: fig_TE_transcribed_vs_mobilized_age.png")
