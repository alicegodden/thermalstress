# =====================================================================
# Active-subset TE analysis (reviewer-requested reframe) - PRIVATE insertions
#   Instead of the genome-wide Kimura landscape (dominated by ancient,
#   fixed copies), focus on the ACTIVE subset: non-reference TE insertions
#   PRIVATE to a single fish (RetroSeq, called per fish; absent from every
#   other fish of the same tissue within 100 bp), summed over the 6 fish
#   of each group.
#
#   RetroSeq reports the broad TE CLASS of each insertion, so the analysis
#   is done at class level:
#   (1) class composition of private insertions vs genome-wide background,
#       at 34°C (heat) and 28°C (control), per organ
#   (2) heat vs control class composition (does heat change WHICH classes move?)
#   (3) a repeat landscape reweighted to the heat-exposed private set's classes
#
#   Class = first element of the RetroSeq MEINFO name (e.g. "LTR-DNA-hybrid"
#   counts as LTR), matching the Python pipeline.
#
# Author: Dr Alice M. Godden
# =====================================================================
library(dplyr); library(tidyr); library(ggplot2); library(reshape2)

CLASSES   <- c("DNA", "LINE", "LTR", "SINE", "RC", "Satellite")
CLASS_COL <- c(DNA = "#4c1d4b", LINE = "#812581", LTR = "#b63679",
               SINE = "#e75263", RC = "#f8850f", Satellite = "#fac228")

# ---------------------------------------------------------------------
# INPUTS
#   Private insertions: either read the private-only VCFs written by
#   private_TE_per_fish_interim.py (PRIVATE_ROOT), or fall back to the
#   summed counts below (from private_TE_by_treatment.tsv).
# ---------------------------------------------------------------------
PRIVATE_ROOT <- "/Users/alicegodden/PycharmProjects/heatstressvariants/miranda/vcf/private_uniq"

private_counts_fallback <- list(
  Ovary  = list(`34C` = c(DNA = 893,  LINE = 96,  LTR = 193, SINE = 53, RC = 21, Satellite = 15),
                `28C` = c(DNA = 858,  LINE = 103, LTR = 206, SINE = 46, RC = 13, Satellite = 26)),
  Testis = list(`34C` = c(DNA = 1118, LINE = 132, LTR = 244, SINE = 89, RC = 16, Satellite = 39),
                `28C` = c(DNA = 756,  LINE = 59,  LTR = 149, SINE = 50, RC = 9,  Satellite = 23)))

divsum_files <- c(Ovary = "FT_all.divsum.csv", Testis = "MT_all.divsum.csv")
private_dirs <- list(Ovary  = c(`34C` = "ovaries/treatment", `28C` = "ovaries/control"),
                     Testis = c(`34C` = "testes/treatment",  `28C` = "testes/control"))

# ---------------------------------------------------------------------
# 1. PRIVATE INSERTIONS - class counts from per-fish private VCFs
# ---------------------------------------------------------------------
te_class <- function(name) {
  first <- toupper(strsplit(name, "-", fixed = TRUE)[[1]][1])
  switch(first, DNA = "DNA", LINE = "LINE", LTR = "LTR", SINE = "SINE",
         RC = "RC", SATELLITE = "Satellite", NA_character_)
}
vcf_class_counts <- function(paths) {
  cnt <- setNames(integer(length(CLASSES)), CLASSES)
  for (p in paths) {
    v <- readLines(p); v <- v[!grepl("^#", v)]
    for (line in v) {
      f <- strsplit(line, "\t")[[1]]
      if (!(f[1] %in% as.character(1:25))) next
      mei <- sub(".*MEINFO=([^,;]+).*", "\\1", f[8])
      cl  <- te_class(mei)
      if (!is.na(cl)) cnt[cl] <- cnt[cl] + 1L
    }
  }
  cnt
}
get_private_counts <- function(organ, temp) {
  d <- file.path(PRIVATE_ROOT, private_dirs[[organ]][[temp]])
  files <- if (dir.exists(d)) list.files(d, pattern = "\\.private\\.vcf$", full.names = TRUE) else character(0)
  if (length(files)) {
    message(sprintf("%s %s: reading %d private VCFs", organ, temp, length(files)))
    return(vcf_class_counts(files))
  }
  message(sprintf("%s %s: private VCFs not found - using summed counts in script", organ, temp))
  private_counts_fallback[[organ]][[temp]][CLASSES]
}

# ---------------------------------------------------------------------
# 2. BACKGROUND - genome-wide bp per class, and full landscape, from divsum
# ---------------------------------------------------------------------
read_divsum <- function(path) {
  d <- read.csv(path, sep = " ", check.names = FALSE)
  d[, !grepl("^Unnamed|ARTEFACT", colnames(d)), drop = FALSE]
}
col_to_class <- function(cols) sapply(cols, function(c) {
  if (c == "Div") return(NA)
  for (cl in CLASSES) if (toupper(c) == toupper(cl) ||
                          startsWith(toupper(c), paste0(toupper(cl), "/"))) return(cl)
  NA
})
background_bp <- function(divsum) {
  cls <- col_to_class(colnames(divsum))
  sapply(CLASSES, function(cl) sum(as.matrix(divsum[, which(cls == cl), drop = FALSE])))
}

# =====================================================================
# RUN - per organ
# =====================================================================
enrich_all <- list(); landscape_all <- list(); heat_vs_ctrl <- list()

for (organ in names(divsum_files)) {
  divsum <- read_divsum(divsum_files[[organ]])
  bg  <- background_bp(divsum); bgp <- bg / sum(bg)
  cnt <- list(`34C` = get_private_counts(organ, "34C"), `28C` = get_private_counts(organ, "28C"))
  
  # (1) private composition vs genome-wide background, each temperature
  for (temp in names(cnt)) {
    a  <- cnt[[temp]]; ap <- a / sum(a)
    chi <- suppressWarnings(chisq.test(x = as.integer(a), p = bgp))
    enrich_all[[paste(organ, temp)]] <- data.frame(
      Organ = organ, Temperature = sub("C$", "°C", temp), Class = CLASSES,
      private_n = as.integer(a), private_pct = 100 * as.numeric(ap),
      genome_pct = 100 * as.numeric(bgp), enrichment = as.numeric(ap / bgp),
      chisq_vs_genome_p = chi$p.value)
  }
  
  # (2) heat vs control class composition (classes with < 5 expected merged into "Other")
  tab <- rbind(`34C` = cnt[["34C"]], `28C` = cnt[["28C"]])
  small <- colSums(tab) < 20
  if (any(small)) tab <- cbind(tab[, !small, drop = FALSE], Other = rowSums(tab[, small, drop = FALSE]))
  ct <- suppressWarnings(chisq.test(tab))
  fi <- tryCatch(fisher.test(tab, simulate.p.value = TRUE, B = 1e5)$p.value, error = function(e) NA)
  heat_vs_ctrl[[organ]] <- data.frame(Organ = organ, n_34C = sum(cnt[["34C"]]), n_28C = sum(cnt[["28C"]]),
                                      classes_tested = paste(colnames(tab), collapse = ","),
                                      chisq = round(unname(ct$statistic), 2), df = unname(ct$parameter),
                                      chisq_p = signif(ct$p.value, 3), fisher_sim_p = signif(fi, 3))
  
  # (3) landscape reweighted by the heat-exposed private set's class enrichment
  ap34 <- cnt[["34C"]] / sum(cnt[["34C"]])
  cls <- col_to_class(colnames(divsum))
  kd  <- melt(divsum, id = "Div", variable.name = "subfam", value.name = "bp")
  kd$Class <- cls[as.character(kd$subfam)]
  kd  <- kd[!is.na(kd$Class), ]
  landscape_all[[organ]] <- kd %>% group_by(Div, Class) %>%
    summarise(bp = sum(bp), .groups = "drop") %>%
    left_join(data.frame(Class = CLASSES, w = as.numeric(ap34 / bgp)), by = "Class") %>%
    mutate(active_weighted = bp * w, Organ = organ)
}

enrich_tab <- bind_rows(enrich_all)
hvc_tab    <- bind_rows(heat_vs_ctrl)
cat("\n===== Private insertions vs genome-wide background (class composition) =====\n")
print(enrich_tab %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame())
cat("\n===== Heat (34°C) vs control (28°C) class composition of private insertions =====\n")
print(hvc_tab)
write.csv(enrich_tab, "private_TE_class_enrichment.csv", row.names = FALSE)
write.csv(hvc_tab,    "private_TE_class_heat_vs_control.csv", row.names = FALSE)

theme_fig <- theme_classic() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text  = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"))

# =====================================================================
# FIGURE 1 - class enrichment of private insertions vs genome (log2 axis)
# =====================================================================
enrich_tab$Class <- factor(enrich_tab$Class, levels = CLASSES)
enrich_tab$Temperature <- factor(enrich_tab$Temperature, levels = c("28°C", "34°C"))
hvc_lab <- hvc_tab %>% mutate(label = sprintf("34 vs 28°C composition: χ² P = %s", signif(chisq_p, 2)))

p1 <- ggplot(enrich_tab, aes(Class, enrichment, fill = Temperature)) +
  geom_col(position = position_dodge(.7), width = .6) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_text(data = hvc_lab, aes(x = 0.6, y = Inf, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 1.5, size = 3.3, fontface = "bold") +
  facet_wrap(~ Organ) +
  scale_fill_manual(values = c(`28°C` = "#2E6DA4", `34°C` = "#C0392B")) +
  scale_y_continuous(trans = "log2", breaks = c(.125, .25, .5, 1, 2, 4)) +
  labs(x = "TE class", y = "Enrichment among private insertions\n(private % / genome-wide %)",
       title = "Class composition of private TE insertions relative to the genome",
       subtitle = "Insertions private to a single fish, summed over 6 fish per group") +
  theme_fig
ggsave("fig_private_TE_enrichment.png", p1, width = 11, height = 5.5, dpi = 600)
ggsave("fig_private_TE_enrichment.pdf", p1, width = 11, height = 5.5)
print(p1)

# =====================================================================
# FIGURE 2 - active-subset repeat landscape (class-reweighted Kimura),
# weighted by the 34°C private set
# =====================================================================
land <- bind_rows(landscape_all)
land$Class <- factor(land$Class, levels = CLASSES)
p2 <- ggplot(land, aes(Div, active_weighted, fill = Class)) +
  geom_col(width = 1) +
  facet_wrap(~ Organ, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = CLASS_COL) +
  coord_cartesian(xlim = c(0, 30)) +
  labs(x = "Kimura substitution level (CpG-adjusted)",
       y = "Activity-weighted TE content (bp × class enrichment)",
       title = "Active-subset repeat landscape",
       subtitle = "Genome-wide Kimura landscape reweighted by TE-class enrichment of private insertions at 34°C") +
  theme_fig
ggsave("fig_private_TE_landscape.png", p2, width = 10, height = 8, dpi = 600)
ggsave("fig_private_TE_landscape.pdf", p2, width = 10, height = 8)
print(p2)

message("Done: fig_private_TE_enrichment.png/.pdf, fig_private_TE_landscape.png/.pdf, ",
        "private_TE_class_enrichment.csv, private_TE_class_heat_vs_control.csv")
