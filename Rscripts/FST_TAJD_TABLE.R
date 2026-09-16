# =====================================================================
# Genes overlapping differentiated regions — per-sex chromosome sets.
#   Female (FT vs FC): Chr 18 and Chr 8
#   Male   (MT vs MC): Chr 24 and Chr 5
# BASE-R ONLY: no flextable / officer / any compiled package.
# Writes:
#   selection_regions_summary_4chr.csv
#   Supp_Table_genes_selection_4chr.tsv   (import into Word/Excel)
#   Supp_Table_genes_selection_4chr.html  (open in browser -> copy into Word)
# Requires gene files: chr18_genes.txt, chr8_genes.txt, chr24_genes.txt,
#                      chr5_genes.txt   (cols: chrom start end attribute)
# =====================================================================

WIN <- 50000; MIN_SNPS <- 50
GAP_MERGE <- 100000; GW_TOP <- 0.99; CHR_TOP <- 0.95
CHR_FEM <- c("18","8")     # ovary panels
CHR_MAL <- c("24","5")     # testis panels
ALL_CHR <- unique(c(CHR_FEM, CHR_MAL))
REPRO <- tolower(c(
  "nanos3","nanos1","nanos2","dnd1","dazl","ddx4","vasa","piwil1","piwil2",
  "dmrt1","sox9a","sox9b","amh","cyp19a1a","figla","foxl2","zar1","bmp15",
  "gdf9","tdrd1","tdrd6","tdrd7","henmt1","mov10l1","ziwi","zili",
  "hsp70","hsp70l","hspa8b","hsp90aa1.2","hsp90ab1","hsf1","hspb1",
  "hsbp1","hsbp1b","dnajb1","dnaja1","dmc1","sycp1","sycp3","spo11",
  "mlh1","rad51","tssk6","boule","gsdf","nanos","buc","rbpms2","ca15b",
  "plcg2","cdh13","grm7","grm7p","parp12","srgap1"))

load_taj <- function(path, lab) {
  d <- read.table(path, header=TRUE, sep="\t")
  data.frame(CHROM=as.character(d$CHROM), WIN=d$BIN_START %/% WIN,
             n=d$N_SNPS, D=suppressWarnings(as.numeric(d$TajimaD)))
}
load_fst <- function(path) {
  d <- read.table(path, header=TRUE, sep="\t")
  data.frame(CHROM=as.character(d$CHROM), WIN=d$BIN_START %/% WIN, FST=d$WEIGHTED_FST)
}
build <- function(fst, tT, tC) {
  f <- load_fst(fst); t <- load_taj(tT,"T"); c <- load_taj(tC,"C")
  m <- merge(f, t, by=c("CHROM","WIN")); m <- merge(m, c, by=c("CHROM","WIN"))
  names(m)[names(m)=="n.x"]<-"nT"; names(m)[names(m)=="D.x"]<-"D_T"
  names(m)[names(m)=="n.y"]<-"nC"; names(m)[names(m)=="D.y"]<-"D_C"
  low <- m$nT < MIN_SNPS | m$nC < MIN_SNPS
  m$D_T[low]<-NA; m$D_C[low]<-NA; m$dD <- m$D_T - m$D_C; m
}
load_genes <- function(path) {
  g <- read.table(path, header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE, fill=TRUE)
  names(g)<-trimws(names(g))
  g <- g[!is.na(g$attribute), ]                     # drop NaN gene names first
  g$attribute <- trimws(as.character(g$attribute))
  g$chrom <- gsub("chr","",as.character(g$chrom))
  g[g$attribute!="" & g$attribute!="NA", c("chrom","start","end","attribute")]
}
sig <- function(dD) {
  if (is.na(dD)) return("differentiated (D not estimable)")
  if (dD < -1) return("directional / sweep under heat (D falls)")
  if (dD >  1) return("balancing / standing variation (D rises)")
  "differentiated, weak D shift"
}
merge_regions <- function(hi) {
  if (nrow(hi)==0) return(list())
  hi <- hi[order(hi$WIN),]; hi$s <- hi$WIN*WIN; hi$e <- hi$WIN*WIN + WIN
  regs<-list(); cur<-NULL
  for (i in seq_len(nrow(hi))) {
    r<-hi[i,]
    if (is.null(cur)) cur<-list(s=r$s,e=r$e,rows=r)
    else if (r$s - cur$e <= GAP_MERGE) { cur$e<-r$e; cur$rows<-rbind(cur$rows,r) }
    else { regs[[length(regs)+1]]<-cur; cur<-list(s=r$s,e=r$e,rows=r) }
  }
  regs[[length(regs)+1]]<-cur; regs
}
annotate <- function(m, sexlab, genes, chrom_set) {
  gw <- quantile(m$FST, GW_TOP, na.rm=TRUE); out<-list()
  for (ch in chrom_set) {
    sub <- m[m$CHROM==ch,]; if (!nrow(sub)) next
    if (sum(sub$FST > gw, na.rm=TRUE)==0) {
      thr <- quantile(sub$FST, CHR_TOP, na.rm=TRUE); otype <- "per-chromosome (top 5%)"
    } else { thr <- gw; otype <- "genome-wide (top 1%)" }
    hi <- sub[sub$FST > thr & !is.na(sub$FST),]
    gch <- genes[genes$chrom==ch,]
    for (reg in merge_regions(hi)) {
      rr<-reg$rows; peak<-rr[which.max(rr$FST),]
      gsub <- gch[gch$start <= reg$e & gch$end >= reg$s,]
      gsub <- gsub[!duplicated(gsub$attribute),]
      if (!nrow(gsub)) next
      cand <- gsub$attribute[tolower(gsub$attribute) %in% REPRO]
      out[[length(out)+1]] <- data.frame(
        Sex=sexlab, Chr=ch,
        Region=sprintf("%.2f-%.2f Mb", reg$s/1e6, reg$e/1e6),
        Outlier_type=otype,
        Peak_FST=round(peak$FST,3),
        TajD_control=round(peak$D_C,2),
        TajD_temperature=round(peak$D_T,2),
        dD_peak=round(peak$dD,2),
        Selection=sig(peak$dD),
        N_genes=nrow(gsub),
        Candidate_genes=if (length(cand)) paste(cand,collapse=", ") else "-",
        All_genes=paste(gsub$attribute,collapse="; "),
        stringsAsFactors=FALSE)
    }
  }
  do.call(rbind, out)
}

genes <- do.call(rbind, lapply(ALL_CHR, function(c) load_genes(sprintf("chr%s_genes.txt", c))))
mal <- build("fst_outputMTvMC.windowed.weir.fst","Taj_2026_MT.Tajima.D","Taj_2026_MC.Tajima.D")
fem <- build("fst_outputFTvFC.windowed.weir.fst","Taj_2026_FT.Tajima.D","Taj_2026_FC.Tajima.D")
summ <- rbind(annotate(fem,"Female (FT vs FC)",genes,CHR_FEM),
              annotate(mal,"Male (MT vs MC)",genes,CHR_MAL))
summ <- summ[order(summ$Sex, summ$Chr, -summ$Peak_FST),]
write.csv(summ, "selection_regions_summary_4chr.csv", row.names=FALSE)
write.table(summ, "Supp_Table_genes_selection_4chr.tsv", sep="\t", row.names=FALSE, quote=FALSE)

# ---------- HTML ----------
esc <- function(x){ x<-gsub("&","&amp;",x); x<-gsub("<","&lt;",x); gsub(">","&gt;",x) }
gene_cell <- function(cand, all) {
  gl <- strsplit(all, "; ")[[1]]
  if (cand != "-") {
    cg <- strsplit(cand, ", ")[[1]]; rest <- setdiff(gl, cg)
    paste0("<b style='color:#a5502a'>", esc(cand), "</b>",
           if (length(rest)) paste0(" &middot; ", esc(paste(rest, collapse=", "))) else "")
  } else esc(paste(gl, collapse=", "))
}
sig_col <- function(s) if (grepl("sweep|directional", s)) "#a5502a" else
  if (grepl("balancing|standing", s)) "#2c5f9e" else "#666666"

rows_html <- character(nrow(summ)); last_sex<-""; shade<-FALSE
for (i in seq_len(nrow(summ))) {
  r <- summ[i,]
  if (r$Sex != last_sex) { shade<-FALSE; last_sex<-r$Sex }
  is_cand <- r$Candidate_genes != "-"
  bg <- if (is_cand) "#fbf0da" else if (shade) "#eef2f7" else "#ffffff"
  dd <- if (is.na(r$dD_peak)) "&ndash;" else sprintf("%+.2f", r$dD_peak)
  dc <- if (is.na(r$TajD_control)) "&ndash;" else sprintf("%.2f", r$TajD_control)
  dt <- if (is.na(r$TajD_temperature)) "&ndash;" else sprintf("%.2f", r$TajD_temperature)
  gw <- grepl("genome-wide", r$Outlier_type)
  ot <- if (gw) "&#9733; GW" else "&#9734; chr"
  td <- function(x, extra="") sprintf("<td style='padding:3px 6px;border:1px solid #ddd;%s'>%s</td>", extra, x)
  rows_html[i] <- paste0("<tr style='background:", bg, "'>",
    td(esc(r$Sex), "font-weight:bold;color:#1f3b63"),
    td(r$Chr, "text-align:center"),
    td(r$Region, "text-align:center"),
    td(ot, "text-align:center;font-size:10px;color:#777"),
    td(sprintf("%.3f", r$Peak_FST), paste0("text-align:center", if (r$Peak_FST>=0.25) ";font-weight:bold" else "")),
    td(dc, "text-align:center;color:#2c5f9e"),
    td(dt, "text-align:center;color:#c85a28"),
    td(dd, paste0("text-align:center", if (!is.na(r$dD_peak) && abs(r$dD_peak)>=1) ";font-weight:bold" else "")),
    td(esc(r$Selection), paste0("color:", sig_col(r$Selection))),
    td(gene_cell(r$Candidate_genes, r$All_genes)), "</tr>")
  shade <- !shade
}
hdr <- c("Sex (contrast)","Chr","Region","Outlier","Peak F<sub>ST</sub>","Taj D control",
         "Taj D 34&deg;C","&Delta;D","Selection signature",
         "Candidate genes (repro/stress bold) &middot; all genes")
hdr_html <- paste0("<th style='padding:4px 6px;border:1px solid #bbb;background:#2c5f9e;color:#fff;'>",
                   hdr, "</th>", collapse="")
caption <- paste0(
  "<b>Supplementary Table S#. Genes overlapping differentiated genomic regions.</b> ",
  "Outlier windows of genetic differentiation (weighted F<sub>ST</sub>, temperature vs control) were identified ",
  "for each sex on chromosomes 18 and 8 (females, FT vs FC) and 24 and 5 (males, MT vs MC), and merged into regions. ",
  "For each region we report the peak-F<sub>ST</sub> window, whether it is a genome-wide (&#9733; GW, top 1%) or ",
  "per-chromosome (&#9734; chr, top 5%) outlier, the paired within-group Tajima's D for control (28&deg;C) and ",
  "temperature (34&deg;C) fish, their difference (&Delta;D = D<sub>34&deg;C</sub> &minus; D<sub>28&deg;C</sub>), the ",
  "inferred selection signature, and overlapping genes. Reproductive/germline/heat-stress candidates are in bold. ",
  "Genes are positional candidates within differentiated regions, not demonstrated targets of selection. Windows: ",
  "50 kb. F<sub>ST</sub>: Weir &amp; Cockerham (1984). D: Tajima (1989).")
html <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'>",
  "<style>body{font-family:Calibri,Arial,sans-serif;font-size:11px;margin:20px}",
  "table{border-collapse:collapse;width:100%}",
  "p.cap{font-size:11px;color:#333;max-width:1150px;line-height:1.4}</style></head><body>",
  "<p class='cap'>", caption, "</p>",
  "<table><thead><tr>", hdr_html, "</tr></thead><tbody>",
  paste(rows_html, collapse=""), "</tbody></table>",
  "<p class='cap'><b>Interpretation.</b> Negative &Delta;D (D falls under heat) is consistent with a directional ",
  "selective sweep in the temperature group; positive &Delta;D (D rises) with balancing selection / standing ",
  "variation; |&Delta;D| &lt; 1 with differentiation lacking a clear diversity shift.</p></body></html>")
writeLines(html, "Supp_Table_genes_selection_4chr.html")
message("Wrote selection_regions_summary_4chr.csv, .tsv and .html (", nrow(summ), " regions)")
