# =====================================================================
# Supplementary table: genes overlapping differentiated regions
#   Chr 18 & Chr 22, both sexes.
#   FST (Temp vs Control, between-group) + Tajima's D (within-group)
#   -> merged outlier regions -> gene overlap -> selection signature
# Author: Dr Alice M. Godden
#
# Outputs:
#   genes_under_selection_chr18_chr22.csv   (full per-gene)
#   selection_regions_summary.csv           (one row per region)
#   Supp_Table_genes_selection.docx         (formatted supplementary table)
#
# Requires: dplyr, tidyr, purrr, flextable, officer
#   install.packages(c("dplyr","tidyr","purrr","flextable","officer"))
# =====================================================================
library(dplyr)
library(tidyr)
library(purrr)
library(flextable)
library(officer)

# ---------------- parameters ----------------------------------------
WIN        <- 50000
CHROMS     <- c("18", "22")
MIN_SNPS   <- 50        # windows below this (either group) unreliable for D
GAP_MERGE  <- 100000    # merge outlier windows within this gap into one region
GW_TOP     <- 0.99      # genome-wide FST outlier quantile
CHR_TOP    <- 0.95      # per-chromosome fallback when a cell has no GW outliers

# curated reproductive / germline / heat-stress candidate genes
REPRO <- tolower(c(
  "nanos3","nanos1","nanos2","dnd1","dazl","ddx4","vasa","piwil1","piwil2",
  "dmrt1","sox9a","sox9b","amh","cyp19a1a","figla","foxl2","zar1","bmp15",
  "gdf9","tdrd1","tdrd6","tdrd7","henmt1","mov10l1","ziwi","zili",
  "hsp70","hsp70l","hspa8b","hsp90aa1.2","hsp90ab1","hsf1","hspb1",
  "hsbp1","hsbp1b","dnajb1","dnaja1",
  "dmc1","sycp1","sycp3","spo11","mlh1","rad51","tssk6","boule",
  "gsdf","nanos","buc","rbpms2","ca15b",
  "plcg2","cdh13","grm7","grm7p","parp12","srgap1"))

# ---------------- loaders -------------------------------------------
load_taj <- function(path, lab) {
  read.table(path, header = TRUE, sep = "\t") %>%
    transmute(CHROM = as.character(CHROM),
              WIN   = BIN_START %/% WIN,
              n     = N_SNPS,
              D     = suppressWarnings(as.numeric(TajimaD))) %>%
    rename(!!paste0("n_", lab) := n, !!paste0("D_", lab) := D)
}
load_fst <- function(path) {
  read.table(path, header = TRUE, sep = "\t") %>%
    transmute(CHROM = as.character(CHROM),
              WIN   = BIN_START %/% WIN,
              FST   = WEIGHTED_FST)
}
build <- function(fst_file, tajT, tajC) {
  m <- load_fst(fst_file) %>%
    inner_join(load_taj(tajT, "T"), by = c("CHROM", "WIN")) %>%
    inner_join(load_taj(tajC, "C"), by = c("CHROM", "WIN"))
  lowp <- m$n_T < MIN_SNPS | m$n_C < MIN_SNPS
  m$D_T[lowp] <- NA; m$D_C[lowp] <- NA
  mutate(m, dD = D_T - D_C)
}
load_genes <- function(path) {
  g <- read.table(path, header = TRUE, sep = "\t", quote = "",
                  stringsAsFactors = FALSE, fill = TRUE)
  names(g) <- trimws(names(g))
  g$attribute <- trimws(as.character(g$attribute))
  g$chrom <- gsub("chr", "", as.character(g$chrom))
  g %>% filter(!is.na(attribute), attribute != "", attribute != "NA") %>%
    select(chrom, start, end, attribute)
}

# ---------------- selection-signature label -------------------------
signature <- function(dD) {
  if (is.na(dD)) return("differentiated (D not estimable)")
  if (dD < -1)  return("directional / sweep under heat (D falls)")
  if (dD >  1)  return("balancing / standing variation (D rises)")
  "differentiated, weak D shift"
}

# ---------------- merge outlier windows into regions ----------------
merge_regions <- function(hi) {
  if (nrow(hi) == 0) return(list())
  hi <- arrange(hi, WIN)
  hi$start <- hi$WIN * WIN; hi$end <- hi$WIN * WIN + WIN
  regions <- list(); cur <- NULL
  for (i in seq_len(nrow(hi))) {
    r <- hi[i, ]
    if (is.null(cur)) {
      cur <- list(start = r$start, end = r$end, rows = r)
    } else if (r$start - cur$end <= GAP_MERGE) {
      cur$end <- r$end; cur$rows <- bind_rows(cur$rows, r)
    } else {
      regions[[length(regions) + 1]] <- cur
      cur <- list(start = r$start, end = r$end, rows = r)
    }
  }
  regions[[length(regions) + 1]] <- cur
  regions
}

# ---------------- annotate one sex ----------------------------------
annotate_sex <- function(m, sexlab, genes) {
  gw99 <- quantile(m$FST, GW_TOP, na.rm = TRUE)
  out <- list()
  for (ch in CHROMS) {
    sub <- filter(m, CHROM == ch)
    if (nrow(sub) == 0) next
    thr <- if (sum(sub$FST > gw99, na.rm = TRUE) == 0)
      min(gw99, quantile(sub$FST, CHR_TOP, na.rm = TRUE)) else gw99
    hi <- filter(sub, FST > thr)
    gch <- filter(genes, chrom == ch)
    for (reg in merge_regions(hi)) {
      rr   <- reg$rows
      peak <- rr[which.max(rr$FST), ]
      pdD  <- peak$dD
      gsub <- gch %>% filter(start <= reg$end, end >= reg$start) %>%
        distinct(attribute, .keep_all = TRUE)
      if (nrow(gsub) == 0) next
      out[[length(out) + 1]] <- tibble(
        Sex = sexlab, Chr = ch,
        Region = sprintf("%.2f-%.2f Mb", reg$start/1e6, reg$end/1e6),
        N_windows = nrow(rr), Genes_in_region = nrow(gsub),
        Gene = gsub$attribute, Gene_start = gsub$start, Gene_end = gsub$end,
        Peak_FST = round(peak$FST, 3),
        TajD_control = round(peak$D_C, 2),
        TajD_temperature = round(peak$D_T, 2),
        dD_peak = round(pdD, 2),
        dD_region_mean = round(mean(rr$dD, na.rm = TRUE), 2),
        Selection = signature(pdD),
        Reproductive_stress_candidate =
          ifelse(tolower(gsub$attribute) %in% REPRO, "yes", ""))
    }
  }
  bind_rows(out)
}

# ---------------- run ------------------------------------------------
genes <- bind_rows(lapply(CHROMS, function(c) load_genes(sprintf("chr%s_genes.txt", c))))

mal <- build("fst_outputMTvMC_windowed_weir.fst",
             "Taj_2026_MT_Tajima.D", "Taj_2026_MC_Tajima.D")
fem <- build("fst_outputFTvFC_windowed_weir.fst",
             "Taj_2026_FT_Tajima.D", "Taj_2026_FC_Tajima.D")

full <- bind_rows(annotate_sex(fem, "Female (FT vs FC)", genes),
                  annotate_sex(mal, "Male (MT vs MC)",   genes))
write.csv(full, "genes_under_selection_chr18_chr22.csv", row.names = FALSE)

# ---------------- region-level summary ------------------------------
summ <- full %>%
  group_by(Sex, Chr, Region) %>%
  summarise(
    Peak_FST = first(Peak_FST),
    TajD_control = first(TajD_control),
    TajD_temperature = first(TajD_temperature),
    dD_peak = first(dD_peak),
    Selection = first(Selection),
    N_genes = n_distinct(Gene),
    Candidate_genes = {
      cg <- unique(Gene[Reproductive_stress_candidate == "yes"])
      if (length(cg)) paste(cg, collapse = ", ") else "-"
    },
    All_genes = paste(unique(Gene), collapse = "; "),
    .groups = "drop") %>%
  arrange(Sex, Chr, desc(Peak_FST))
write.csv(summ, "selection_regions_summary.csv", row.names = FALSE)

# ---------------- formatted Word table (flextable) ------------------
# build a display gene column: candidates first, then the rest
disp <- summ %>%
  mutate(Genes = mapply(function(cand, all) {
    genes_all <- strsplit(all, "; ")[[1]]
    if (cand != "-") {
      rest <- setdiff(genes_all, strsplit(cand, ", ")[[1]])
      paste0(cand, if (length(rest)) paste0("  ·  ", paste(rest, collapse = ", ")) else "")
    } else paste(genes_all, collapse = ", ")
  }, Candidate_genes, All_genes),
  dD = sprintf("%+.2f", dD_peak)) %>%
  select(Sex, Chr, Region, Peak_FST, TajD_control, TajD_temperature,
         dD, Selection, Genes)

# flag: TRUE where a region contains a curated candidate gene (used for shading)
cand_flag <- summ$Candidate_genes != "-"

ft <- flextable(disp) %>%
  set_header_labels(
    Sex = "Sex (contrast)", Chr = "Chr", Region = "Region",
    Peak_FST = "Peak FST", TajD_control = "Taj D control",
    TajD_temperature = "Taj D 34 \u00b0C", dD = "\u0394D",
    Selection = "Selection signature",
    Genes = "Candidate genes (repro/stress bold) \u00b7 all genes") %>%
  # colour the two D columns
  color(j = "TajD_control", color = "#2c5f9e", part = "body") %>%
  color(j = "TajD_temperature", color = "#c85a28", part = "body") %>%
  bold(j = "Sex", part = "body") %>%
  color(j = "Sex", color = "#1f3b63", part = "body") %>%
  # header styling
  bg(bg = "#2c5f9e", part = "header") %>%
  color(color = "white", part = "header") %>%
  bold(part = "header") %>%
  # shade the signature cell of rows whose region contains a candidate gene
  bg(i = cand_flag, j = "Selection", bg = "#fbf0da") %>%
  fontsize(size = 8, part = "all") %>%
  align(j = c("Chr","Region","Peak_FST","TajD_control",
              "TajD_temperature","dD"), align = "center", part = "all") %>%
  border_outer(border = fp_border(color = "#bbbbbb", width = 1)) %>%
  border_inner_h(border = fp_border(color = "#dddddd", width = .5)) %>%
  width(j = "Sex", width = 1.3) %>%
  width(j = "Chr", width = 0.4) %>%
  width(j = "Region", width = 1.1) %>%
  width(j = c("Peak_FST","TajD_control","TajD_temperature","dD"), width = 0.7) %>%
  width(j = "Selection", width = 1.9) %>%
  width(j = "Genes", width = 4.4)

# helper: none needed — cand_flag computed above from summ

caption <- paste0(
  "Supplementary Table S#. Genes overlapping differentiated genomic regions on ",
  "chromosomes 18 and 22. Outlier windows of genetic differentiation (weighted FST, ",
  "temperature vs control) were identified for each sex and merged into regions. For each ",
  "region we report the peak-FST window, the paired within-group Tajima's D at that window ",
  "for control (28 \u00b0C) and temperature (34 \u00b0C) fish, their difference ",
  "(\u0394D = D34\u00b0C - D28\u00b0C), the inferred selection signature, and overlapping genes. ",
  "Reproductive/germline/heat-stress candidates are in bold. Genes are positional candidates ",
  "within differentiated regions, not demonstrated targets of selection. Windows: 50 kb. ",
  "FST: Weir & Cockerham (1984). D: Tajima (1989).")

doc <- read_docx() %>%
  body_add_par(caption, style = "Normal") %>%
  body_add_par("") %>%
  body_add_flextable(ft) %>%
  body_end_section_landscape()          # landscape page for the wide table
print(doc, target = "Supp_Table_genes_selection.docx")

message("Done: genes_under_selection_chr18_chr22.csv, selection_regions_summary.csv, Supp_Table_genes_selection.docx")
