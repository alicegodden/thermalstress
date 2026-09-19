# Rscripts — *thermalstress*

R scripts for the analysis and figures in **Godden *et al.*, "Sex differences in sRNA-mediated transposon regulation under thermal stress in zebrafish germ cells."**

These scripts take the processed outputs of the upstream pipelines (DESeq2 tables, RetroSeq VCFs, RepeatMasker divergence summaries, VCFtools population-genetics output, Telescope count matrices) and produce the statistical analyses and figures in the manuscript. Each script is self-contained and can be run independently once its input files are present in the working directory.

---

## Overview

Adult zebrafish were held at control (28 °C) or elevated (34 °C) temperature for two weeks, and ovaries and testes were profiled by whole-genome DNA-seq, mRNA-seq and small RNA-seq. These scripts cover the transposable element (TE), population-genetics, expression and transgenerational analyses.

---

## Scripts

| Script | Purpose |
|--------|---------|
| `sex_selection_figures.R` | Genome-wide F<sub>ST</sub> (between-group) and Tajima's D (within-group) in 50 kb windows, per sex; genome-wide Manhattan + ΔD contrast, and per-chromosome dual-axis panels with centromere annotation. |
| `TE_counts_revised.R` | De novo (treatment-unique) TE insertion counts per family; per-family GLMMs (Poisson / negative binomial, Temperature vs Control within each sex) and testes-vs-ovaries figures. |
| `telescope_deseq2.R` | DESeq2 on Telescope locus-level TE count matrices (Temperature vs Control), per sex, producing differential TE expression results. |
| `active_te_landscape.R` | Active-subset TE analysis: class composition of temperature-unique (RetroSeq) insertions vs genome-wide background (χ² enrichment), and an activity-reweighted Kimura repeat landscape. |
| `te_open_chromatin_enrichment.R` | Enrichment of temperature-unique de novo insertions in open chromatin (DANIO-CODE DOPES regions) vs genome-wide coverage; binomial test, per sex. |
| `te_open_chromatin_temp_ctrl.R` | As above, including control insertions: open/closed counts for all four groups, enrichment tests, and Temperature-vs-Control comparison (Wilcoxon). |
| `te_chromatin_boxplot.R` | Boxplot of de novo insertions in open vs closed chromatin, ovary vs testis, with statistical tests. |
| `te_transcription_vs_mobilization_age.R` | Compares the evolutionary age (Kimura) of transcriptionally active TEs (Telescope DE) with mobilized TEs (RetroSeq active subset), per sex, with a KS test. |
| `make_table_baseR.R` | Base-R (no external packages) builder for the "genes overlapping differentiated regions" supplementary table (HTML + TSV). |
| `make_table_4chr.R` | Per-sex, multi-chromosome version of the annotated selection table (gene overlap with F<sub>ST</sub>/Tajima's D outlier regions). |

*(Filenames may differ slightly from those in the manuscript's methods; adjust input paths at the top of each script to match your data.)*

---

## Inputs

Depending on the script, you will need some of:

- **VCFtools output** — `*.windowed.weir.fst` (F<sub>ST</sub>) and `*.Tajima.D` (Tajima's D), per group.
- **RetroSeq VCFs** — treatment-unique de novo TE insertion calls (temperature and control, per sex).
- **RepeatMasker divergence summaries** — `*_all.divsum.csv` (Kimura landscape per assembly).
- **Telescope count matrices / DESeq2 tables** — locus-level TE expression.
- **Gene annotation** — `chrN_genes.txt` (columns: `chrom`, `start`, `end`, `attribute`).
- **Open chromatin** — DANIO-CODE DOPES BED (`daniocode_hub_280355_dopes_all.txt`).
- **Centromere positions** — GRCz11, provided inline in the relevant scripts.

Chromosome coordinates are GRCz11; TE windows are 50 kb throughout.

---

## Dependencies

Most scripts use base R plus the tidyverse:

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork",
                   "reshape2", "viridisLite"))
```

Additional packages by script:

- **DE / GLMM:** `lme4`, `performance` (`TE_counts_revised.R`); `DESeq2` (`telescope_deseq2.R`, Bioconductor)
- **Genomic overlap:** `GenomicRanges` (`te_open_chromatin_enrichment.R`, Bioconductor)
- **Volcano plots:** `EnhancedVolcano`, `ggrepel` (Bioconductor / CRAN)

Bioconductor packages:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "GenomicRanges", "EnhancedVolcano"))
```

`make_table_baseR.R` deliberately uses **base R only** (no compiled dependencies) for environments where package installation is restricted.

---

## Usage

Each script is run from the directory containing its input files, e.g.:

```r
source("sex_selection_figures.R")
```

Input filenames are set near the top of each script — edit these to match your files before running. Figures are written to the working directory at publication resolution.

---

## Data availability

Raw sequencing data are archived at the European Nucleotide Archive under accession **PRJEB72689**.

---

## Citation

If you use these scripts, please cite:

> Godden A.M., De Coriolis J.-C., Sullivan A., Drake C., Immler S. *Sex differences in sRNA-mediated transposon regulation under thermal stress in zebrafish germ cells.* (in review).

---

## Contact

Dr. Alice M. Godden — School of Biological Sciences, University of East Anglia
