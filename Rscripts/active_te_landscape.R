# =====================================================================
# Active-subset TE analysis (reviewer-requested reframe)
#   Instead of the genome-wide Kimura landscape (dominated by ancient,
#   fixed copies and near-identical between conditions), we focus on the
#   ACTIVE subset: the de novo insertions unique to the temperature
#   genome (RetroSeq), following the same rationale as the RetroSeq
#   analysis.
#
#   RetroSeq reports the broad TE CLASS of each insertion (not subfamily
#   or per-copy divergence), so the analysis is done at class level:
#   (1) class composition of active insertions vs genome-wide background
#   (2) a repeat landscape reweighted to the active subset's classes.
#
# Author: Dr Alice M. Godden
# =====================================================================
library(dplyr); library(tidyr); library(ggplot2)
library(reshape2); library(viridisLite)

CLASSES <- c("DNA","LINE","LTR","SINE","RC","Satellite")
CLASS_COL <- c(DNA="#4c1d4b", LINE="#812581", LTR="#b63679",
               SINE="#e75263", RC="#f8850f", Satellite="#fac228")

# ---------------------------------------------------------------------
# 1. ACTIVE SET — class counts of temperature-unique insertions (VCF)
# ---------------------------------------------------------------------
vcf_class_counts <- function(vcf_path) {
  v <- readLines(vcf_path); v <- v[!grepl("^#", v)]
  keep <- as.character(1:25)
  cnt <- setNames(integer(length(CLASSES)), CLASSES)
  for (line in v) {
    f <- strsplit(line, "\t")[[1]]
    if (!(f[1] %in% keep)) next
    info <- f[8]
    mei <- sub(".*MEINFO=([^,]+),.*", "\\1", info)   # broad class only
    for (cl in CLASSES) if (toupper(mei) == toupper(cl)) { cnt[cl] <- cnt[cl] + 1; break }
  }
  cnt
}

# ---------------------------------------------------------------------
# 2. BACKGROUND — genome-wide bp per class, and full landscape, from divsum
# ---------------------------------------------------------------------
read_divsum <- function(path) {
  d <- read.csv(path, sep = " ", check.names = FALSE)
  d <- d[, !grepl("^Unnamed|ARTEFACT", colnames(d)), drop = FALSE]
  d
}
# map each divsum column to a broad class
col_to_class <- function(cols) {
  sapply(cols, function(c) {
    if (c == "Div") return(NA)
    hit <- NA
    for (cl in CLASSES) if (toupper(c) == toupper(cl) ||
                            startsWith(toupper(c), paste0(toupper(cl), "/"))) { hit <- cl; break }
    hit
  })
}
background_bp <- function(divsum) {
  cls <- col_to_class(colnames(divsum))
  sapply(CLASSES, function(cl) sum(as.matrix(divsum[, which(cls == cl), drop = FALSE])))
}

# =====================================================================
# RUN — per organ (temperature genome vs its own background)
#   Ovary : FT insertions + FT divsum ; Testis : MT insertions + MT divsum
# =====================================================================
inputs <- list(
  Ovary  = list(vcf = "female_exp_unique_8_win100_filterpy_gq750_fl8.vcf", divsum = "FT_all.divsum.csv"),
  Testis = list(vcf = "male_exp_unique_8_win100_filterpy_gq750_fl8.vcf",   divsum = "MT_all.divsum.csv"))
enrich_all <- list(); landscape_all <- list()

for (organ in names(inputs)) {
  active <- vcf_class_counts(inputs[[organ]]$vcf)
  divsum <- read_divsum(inputs[[organ]]$divsum)
  bg     <- background_bp(divsum)

  # --- enrichment: active % vs background % ---
  ap <- active / sum(active)
  bp <- bg / sum(bg)
  enr <- data.frame(Organ = organ, Class = CLASSES,
                    active_n = as.integer(active),
                    active_pct = 100 * ap, bg_pct = 100 * bp,
                    enrichment = as.numeric(ap / bp))
  # chi-square: active composition vs background
  expct <- bp * sum(active)
  chi <- suppressWarnings(chisq.test(x = as.integer(active), p = bp))
  enr$chisq_p <- chi$p.value
  enrich_all[[organ]] <- enr

  # --- active-subset landscape: reweight genome landscape by class enrichment ---
  cls <- col_to_class(colnames(divsum))
  kd <- melt(divsum, id = "Div", variable.name = "subfam", value.name = "bp")
  kd$Class <- cls[as.character(kd$subfam)]
  kd <- kd[!is.na(kd$Class), ]
  # collapse to class x Div, weight by enrichment (active over-representation)
  wl <- kd %>% group_by(Div, Class) %>% summarise(bp = sum(bp), .groups = "drop") %>%
    left_join(data.frame(Class = CLASSES, w = as.numeric(ap / bp)), by = "Class") %>%
    mutate(active_weighted = bp * w, Organ = organ)
  landscape_all[[organ]] <- wl
}

enrich_tab <- bind_rows(enrich_all)
cat("\n===== Active insertions vs genome-wide background (class composition) =====\n")
print(enrich_tab %>% mutate(across(where(is.numeric), ~round(.x, 3))) %>% as.data.frame())
write.csv(enrich_tab, "active_TE_class_enrichment.csv", row.names = FALSE)

# =====================================================================
# FIGURE 1 — class enrichment among active insertions (log2 fold)
# =====================================================================
p1 <- ggplot(enrich_tab, aes(Class, enrichment, fill = Organ)) +
  geom_col(position = position_dodge(.7), width = .6, alpha = .9) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  scale_fill_manual(values = c(Ovary = "#d98f52", Testis = "#c85a28")) +
  scale_y_continuous(trans = "log2", breaks = c(.25,.5,1,2,3)) +
  labs(x = "TE class", y = "Enrichment among active insertions\n(active % / genome-wide %)",
       title = "Heat-mobilized TEs are compositionally distinct from the genomic background",
       subtitle = "Temperature-unique de novo insertions vs genome-wide TE content") +
  theme_classic() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"),
        axis.title = element_text(size = 13, face = "bold"),
        axis.text  = element_text(size = 12, face = "bold"),
        legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"))
ggsave("fig_active_TE_enrichment.png", p1, width = 10, height = 6, dpi = 600)

p1
# =====================================================================
# FIGURE 2 — active-subset repeat landscape (class-reweighted Kimura)
# genome-wide landscape reweighted by each class's activity, per organ
# =====================================================================
land <- bind_rows(landscape_all)
land$Class <- factor(land$Class, levels = CLASSES)

p2 <- ggplot(land, aes(x = Div, y = active_weighted, fill = Class)) +
  geom_col(width = 1) +
  facet_wrap(~ Organ, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = CLASS_COL) +
  coord_cartesian(xlim = c(0, 30)) +
  labs(x = "Kimura substitution level (CpG-adjusted)",
       y = "Activity-weighted TE content (bp x class enrichment)",
       title = "Active-subset repeat landscape",
       subtitle = "Genome-wide Kimura landscape reweighted by TE-class activity (temperature-unique insertions)") +
  theme_classic() +
  theme(plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text  = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(face = "bold"), legend.title = element_text(face = "bold"))
ggsave("fig_active_TE_landscape.png", p2, width = 10, height = 8, dpi = 600)

message("Done: fig_active_TE_enrichment.png, fig_active_TE_landscape.png, active_TE_class_enrichment.csv")
p2
