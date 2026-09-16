# =====================================================================
# Sex-specific selection signals under thermal stress
#   Between-population : Weighted FST  (Temperature vs Control)
#   Within-population  : Tajima's D    (each group separately)
# Females (FT vs FC) and Males (MT vs MC), 50 kb windows.
# Author: Dr Alice M. Godden
# =====================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)   # for multi-panel assembly

WIN <- 50000
CHROMS <- as.character(1:25)
BLUE <- "#2c5f9e"; ORNG <- "#c85a28"; GOLD <- "#e8b64c"; FSTC <- "#8a8f98"
MIN_SNPS <- 50      # windows below this (in either group) are unreliable for D

# centromere positions (bp), zebrafish GRCz11 — from tajD_fst_genes.R
centromere_data <- data.frame(
  CHROM = as.character(1:25),
  CENTROMERE_POS = c(33133433, 19743539, 49120158, 25975534, 50391699,
                     34866988, 60898697, 30099230, 13283057, 9894335,
                     9507769, 27779269, 19555524, 16933009, 13943685,
                     19463555, 48225166, 24528805, 19777257, 11319545,
                     28908825, 21520944, 13631408, 15196971, 20793603))

# ---- bold theme add-on: makes ALL text bold, axis text/titles included ----
bold_theme <- theme(
  text        = element_text(face = "bold"),
  axis.title  = element_text(face = "bold"),
  axis.text   = element_text(face = "bold"),
  axis.text.x = element_text(face = "bold"),
  axis.text.y = element_text(face = "bold"),
  plot.title  = element_text(face = "bold"),
  legend.text = element_text(face = "bold"),
  legend.title= element_text(face = "bold"),
  strip.text  = element_text(face = "bold")
)

# ---- loaders -------------------------------------------------------
load_taj <- function(path, lab) {
  read.table(path, header = TRUE, sep = "\t") %>%
    transmute(CHROM = as.character(CHROM),
              WIN = BIN_START %/% WIN,
              !!paste0("n_", lab)  := N_SNPS,
              !!paste0("D_", lab)  := suppressWarnings(as.numeric(TajimaD)))
}
load_fst <- function(path) {
  read.table(path, header = TRUE, sep = "\t") %>%
    transmute(CHROM = as.character(CHROM),
              WIN = BIN_START %/% WIN,
              FST = WEIGHTED_FST)
}

build <- function(fst_file, tajT, tajC) {
  m <- load_fst(fst_file) %>%
    inner_join(load_taj(tajT, "T"), by = c("CHROM", "WIN")) %>%
    inner_join(load_taj(tajC, "C"), by = c("CHROM", "WIN")) %>%
    filter(CHROM %in% CHROMS)
  # mask low-power windows for D only
  lowp <- m$n_T < MIN_SNPS | m$n_C < MIN_SNPS
  m$D_T[lowp] <- NA; m$D_C[lowp] <- NA
  m %>% mutate(CHROM = factor(CHROM, levels = CHROMS),
               dD = D_T - D_C,
               pos_mb = WIN * 50000 / 1e6)   # window index x window size (bp) -> Mb
}

fem <- build("fst_outputFTvFC.windowed.weir.fst",
             "Taj_2026_FT.Tajima.D", "Taj_2026_FC.Tajima.D")
mal <- build("fst_outputMTvMC.windowed.weir.fst",
             "Taj_2026_MT.Tajima.D", "Taj_2026_MC.Tajima.D")

# =====================================================================
# FIGURE 1 : genome-wide overview (per sex): FST manhattan + dD contrast
# =====================================================================
genome_x <- function(m) {
  sizes <- m %>% group_by(CHROM) %>% summarise(mx = max(WIN) + 1, .groups = "drop") %>%
    right_join(tibble(CHROM = factor(CHROMS, levels = CHROMS)), by = "CHROM") %>%
    mutate(mx = ifelse(is.na(mx), 0, mx))
  sizes$offset <- cumsum(c(0, head(sizes$mx + 4, -1)))
  m %>% left_join(sizes, by = "CHROM") %>% mutate(gx = offset + WIN)
}

fst_manhattan <- function(m, title) {
  mg <- genome_x(m); thr <- quantile(mg$FST, 0.99, na.rm = TRUE)
  ticks <- genome_x(m) %>% group_by(CHROM) %>% summarise(mid = mean(gx), .groups = "drop")
  mg <- mg %>% mutate(shade = ifelse(as.integer(CHROM) %% 2 == 0, "a", "b"),
                      hi = FST > thr)
  ggplot(mg, aes(gx, FST)) +
    geom_point(aes(color = shade), size = .5, alpha = .55) +
    geom_point(data = filter(mg, hi), color = GOLD, size = 1.8,
               shape = 17, alpha = 0.5) +
    geom_hline(yintercept = thr, linetype = "dashed", color = "#7a5b12", linewidth = .3) +
    scale_color_manual(values = c(a = BLUE, b = "#7d97bd"), guide = "none") +
    scale_x_continuous(breaks = ticks$mid, labels = ticks$CHROM) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(title = title, x = NULL, y = expression(bold("Weighted F"[ST]))) +
    theme_minimal(base_size = 13) +
    bold_theme +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 6, face = "bold"),
          plot.title = element_text(face = "bold", hjust = .5))
}

dD_violin <- function(m) {
  thr99 <- quantile(m$FST, 0.99, na.rm = TRUE)
  thr95 <- quantile(m$FST, 0.95, na.rm = TRUE)
  d <- m %>% filter(!is.na(dD)) %>%
    mutate(grp = case_when(FST > thr99 ~ "high FST (top 1%)",
                           FST <= thr95 ~ "background",
                           TRUE ~ NA_character_)) %>%
    filter(!is.na(grp)) %>%
    mutate(grp = factor(grp, levels = c("background", "high FST (top 1%)")))
  means <- d %>% group_by(grp) %>% summarise(m = mean(dD), .groups = "drop")
  
  # Mann-Whitney U (Wilcoxon rank-sum): is signed dD shifted at high-FST loci?
  bg <- d$dD[d$grp == "background"]; hi <- d$dD[d$grp == "high FST (top 1%)"]
  wt <- wilcox.test(hi, bg)
  r_rb <- as.numeric(1 - 2 * wt$statistic / (length(hi) * length(bg)))  # rank-biserial
  p_txt <- if (wt$p.value < 2.2e-16) "P < 2.2e-16" else sprintf("P = %.2g", wt$p.value)
  lab <- sprintf("Wilcoxon %s\nr = %+.2f,  n = %d", p_txt, r_rb, length(hi))
  
  ggplot(d, aes(grp, dD, fill = grp)) +
    geom_violin(trim = TRUE, alpha = .8, color = "gray50") +
    stat_summary(fun = mean, geom = "crossbar", width = .5, linewidth = .3) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray") +
    geom_text(data = means, aes(grp, 3.4, label = sprintf("%+.2f", m)),
              inherit.aes = FALSE, fontface = "bold") +
    annotate("segment", x = 1, xend = 2, y = 4.15, yend = 4.15, linewidth = .4) +
    annotate("text", x = 1.5, y = 4.55, label = lab, size = 3.1, lineheight = .9,
             fontface = "bold") +
    scale_fill_manual(values = c("background" = "#b9c4d6",
                                 "high FST (top 1%)" = GOLD), guide = "none") +
    coord_cartesian(ylim = c(-4, 5), clip = "off") +
    labs(x = NULL, y = expression(bold(Delta*"D = D"[34*degree*C] - "D"[28*degree*C]))) +
    theme_minimal(base_size = 13) +
    bold_theme +
    theme(panel.grid.minor = element_blank(),
          plot.margin = margin(t = 18, r = 8, b = 6, l = 6))
}

fig1 <- (fst_manhattan(fem, "Females (ovary): 34 vs 28 \u00b0C") |
           fst_manhattan(mal, "Males (testis): 34 vs 28 \u00b0C")) /
  (dD_violin(fem) | dD_violin(mal)) +
  plot_annotation(
    title = "Sex-specific selection signals: between-group FST and within-group Tajima's D",
    subtitle = "FST locates Control-vs-Temperature differentiation; \u0394D shows whether Tajima's D shifts at those loci (negative = sweep under heat). 50 kb windows.",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 11, color = "grey40", face = "bold")))
ggsave("fig_selection_genomewide.png", fig1, width = 16, height = 8, dpi = 600)
ggsave("fig_selection_genomewide.pdf", fig1, width = 16, height = 8)

# =====================================================================
# FIGURE 2 : per-chromosome dual axis, BOTH D tracks + FST
# =====================================================================
FST_SLOPE <- 0.06   # maps FST onto the Tajima's D axis (0-aligned)

perchrom <- function(m, chrom, title) {
  sub <- m %>% filter(CHROM == chrom) %>% arrange(WIN)
  thr99 <- quantile(m$FST, 0.99, na.rm = TRUE)
  sub <- sub %>% mutate(FST_disp = FST / FST_SLOPE, hi = FST > thr99)
  cen_mb <- centromere_data$CENTROMERE_POS[centromere_data$CHROM == chrom] / 1e6
  ggplot(sub, aes(pos_mb)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray") +
    geom_vline(xintercept = cen_mb, linetype = "22", color = "gray30",
               linewidth = .7, alpha = .6) +
    annotate("text", x = cen_mb, y = 4.7, label = "centromere",
             size = 2.6, fontface = "bold", color = "gray30", hjust = .5) +
    geom_point(aes(y = FST_disp), color = FSTC, shape = 15, size = 1, alpha = .45) +
    geom_point(data = filter(sub, hi), aes(y = FST_disp),
               color = GOLD, shape = 17, size = 2.6, alpha = 0.5) +
    geom_point(aes(y = D_C), color = BLUE, size = 1.1, alpha = .7) +
    geom_point(aes(y = D_T), color = ORNG, size = 1.1, alpha = .7) +
    scale_y_continuous(name = "Tajima's D", limits = c(-4, 5),
                       sec.axis = sec_axis(~ . * FST_SLOPE,
                                           name = expression(bold("Weighted F"[ST])))) +
    labs(title = title, x = "Position (Mb)") +
    theme_minimal(base_size = 13) +
    bold_theme +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = .5),
          axis.title.y.left = element_text(color = "#20242b", face = "bold"),
          axis.title.y.right = element_text(color = "#7a5b12", face = "bold"),
          axis.text.y.right = element_text(color = "#7a5b12", face = "bold"))
}

# choose chromosomes with strongest signal (edit as needed)
fig2 <- (perchrom(mal, "24", "Testes \u00b7 Chr 24") | perchrom(mal, "5", "Testes \u00b7 Chr 5")) /
  (perchrom(fem, "18", "Ovaries \u00b7 Chr 18") | perchrom(fem, "8", "Ovaries \u00b7 Chr 8")) +
  plot_annotation(
    title = "Per-chromosome FST with paired Tajima's D (blue = 28 \u00b0C, orange = 34 \u00b0C)",
    subtitle = "Separation of the two D tracks at gold FST peaks reveals the selection signature.",
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 11, color = "grey40", face = "bold")))
ggsave("fig_selection_perchrom.png", fig2, width = 16, height = 9, dpi = 600)
ggsave("fig_selection_perchrom.pdf", fig2, width = 16, height = 9)

# =====================================================================
# candidate table: high-FST windows ranked by |dD|, per sex
# =====================================================================
candidates <- function(m, sexlab) {
  thr99 <- quantile(m$FST, 0.99, na.rm = TRUE)
  m %>% filter(FST > thr99, !is.na(dD)) %>%
    mutate(sex = sexlab, WIN_START = WIN * WIN,
           signature = case_when(dD < -1 ~ "sweep under heat (D drops)",
                                 dD >  1 ~ "standing variation (D rises)",
                                 TRUE ~ "differentiated, weak D shift")) %>%
    arrange(desc(abs(dD))) %>%
    select(sex, CHROM, WIN_START, FST, D_C, D_T, dD, n_C, n_T, signature)
}
write.csv(bind_rows(candidates(mal, "male"), candidates(fem, "female")),
          "selection_candidates.csv", row.names = FALSE)

message("Done: fig_selection_genomewide.*, fig_selection_perchrom.*, selection_candidates.csv")
