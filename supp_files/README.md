# 🧬 The Molecular Response to Heat Stress in Zebrafish Gonads

## 📄 Publication Data Repository

This repository contains the supplementary and raw output files for the paper:
**"The molecular response to heat stress in male and female zebrafish gonads from genome to transcriptome and back."**

---

## 🔬 Data Overview: Multi-Omics Analysis

This study integrates multiple layers of genomic response across two sexes (Testes and Ovaries) following two weeks of treatment (28°C vs. 34°C).

The files are organized by molecular component (Gene, TE, miRNA, piRNA, structRNA).

### 🌡️ Environment and Quality Control

| ID | Filepath | Description |
| :---: | :--- | :--- |
| **1** | `water_testing.xlsx` | Raw data detailing **water quality** parameters collected during the 28°C (Control) and 34°C (Treatment) temperature exposure for two weeks. |

### 📈 Core Differential Expression Analysis (DESeq2)

This section contains the results for differential expression of **Genes** and **Transposable Elements (TEs)** between the treated (T) and control (C) groups.

| ID | Filepath | Description |
| :---: | :--- | :--- |
| **2** | `2-te_tocontrol_deseq2.csv` | **Gene Expression (Testes):** Full DESeq2 output comparing gene expression in male gonads (34°C vs. 28°C). |
| **3** | `3-ov_tocontrol_deseq2.csv` | **Gene Expression (Ovaries):** Full DESeq2 output comparing gene expression in female gonads (34°C vs. 28°C). |
| **4** | `4-te_tetrans_deseq2.csv` | **TE Expression (Testes, TEtranscripts):** DESeq2 output for Transposable Elements using the TEtranscripts quantification method. |
| **5** | `5-ov_tetrans_deseq2.csv` | **TE Expression (Ovaries, TEtranscripts):** DESeq2 output for Transposable Elements using the TEtranscripts quantification method. |

### 🦠 Transposable Element Validation (Telescope)

| ID | Filepath | Description |
| :---: | :--- | :--- |
| **10** | `10-telescope_male_raw.csv.zip` | **TE Expression (Testes, Telescope):** Raw quantification output from the Telescope pipeline for male gonads (ZIPPED). |
| **11** | `11-telescope_female_raw.csv.zip` | **TE Expression (Ovaries, Telescope):** Raw quantification output from the Telescope pipeline for female gonads (ZIPPED). |

### 🔑 Small Non-Coding RNA Regulation (miRNA, piRNA, structRNA)

This section details the differential expression of regulatory small RNAs. Files ending in `volcano.tsv`/`.csv` are typically ready for plotting.

| ID | Filepath | Type | Sex |
| :---: | :--- | :--- | :---: |
| **6** | `6-TESTES_DEDUP_HAIRPIN_TVsC_volcano.tsv` | **Hairpin miRNA** | ♂️ Testes |
| **7** | `7-OVARIES_DEDUP_HAIRPIN_TVsC_volcano.tsv` | **Hairpin miRNA** | ♀️ Ovaries |
| **8** | `8-TESTES_DEDUP_MATURE_TVsC_NA_use2.tsv` | **Mature miRNA** | ♂️ Testes |
| **9** | `9-OVARIES_DEDUP_MATURE_TVsC_NA_use.tsv` | **Mature miRNA** | ♀️ Ovaries |
| **12** | `12-male_tesmall_pirna_volcano.csv` | **piRNA** | ♂️ Testes |
| **13** | `13-female_tesmall_pirna_volcano.tsv` | **piRNA** | ♀️ Ovaries |
| **14** | `14-male_structRNA_tesmall_deseq2.csv` | **Structural RNA** | ♂️ Testes |
| **15** | `15-female_structRNA_tesmall_deseq2.csv` | **Structural RNA** | ♀️ Ovaries |
