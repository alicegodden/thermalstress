# 🐟 Environmental Stress Scripts  
**Thermal stress Zebrafish genomics project**  

A collection of analysis scripts for water quality, genomics, TE insertion studies, and more — developed for thermal stress experiments in zebrafish.  


---

## 📂 Contents  

###  bash scripts  
| Script | Description |
|--------|-------------|
| `bash_scripts.txt` | NF-CORE RNA-seq and sm-RNAseq config files and special bash and R scripts for data generation and analysis |
| `metadata_rnaseq.txt` | Sample metadata with organ and sex for RNA-seq data |
| `metadata_smrnaseq.txt` | Sample metadata with organ and sex for RNA-seq data |

---

### 💧 Water Analysis  
| Script | Description |
|--------|-------------|
| `water_testing.py` | Water testing parameters and graph generation |

------

### piRNA ping-pong Analysis  
| Script | Description |
|--------|-------------|
| `ping_pong.py` | Searching for ping-pong signatures |

---

### 📊 R Scripts (`/Rscripts`)  
| Script | Description |
|--------|-------------|
| `AG_RNAseq.Rmd` | R notebook for DESeq2 analyses |
| `manhattan_99.R` | Manhattan plots with 99th percentile plotting |
| `pheatmap_logcpm_z_heatstress.R` | Normalised counts & heatmaps from DESeq2 outputs |
| `tajD_fst_genes.R` | Tajima's D and FST plots with gene tracks at chromosomal level |
| `telescopeCounts.R` | Custom scripts for analysis & assimilation of count data from Telescope pipeline |
| `testing_glm_HS_retroseq.R` | GLM of Retroseq data at family level |
| `zfish_embryo_models.R` | GLM of embryo and mortality data |

---

### 🐍 Python Scripts  

#### 🧬 SNPEFF Analysis  
| Script | Description |
|--------|-------------|
| `filter.py` | Filtering SnpEff and SnpSift outputs |
| `tstvplot.py` | Plotting Ts/Tv ratios from SnpEff |
| `snpeff.py` | Plotting variants from SnpEff and SnpSift |
| `vcfinput_phenogram.py` | Phenogram plot |
| `autobubble_goplot.py` | GO terms plotting (.csv from ShinyGO) with annotated p-values |
| `unique.py` | Filtering non-unique rows from Retroseq outputs to create unique variants file |


#### 🧪 Retroseq Analysis  
| Script | Description |
|--------|-------------|
| `filter_retroseq.py` | Filter Retroseq outputs on GQ & FL for confident TE insertion calls |
| `retroseq.py` + `vcfinput_phenogram.py` | Plot Retroseq outputs & chromosome shapes |
| `retroseq_phenogram_TEcolours.py` | Group TEs by family & plot in colours by family |
| `retroseq_phenogram_TEfam_openChrom.py` | Same as above, with chromatin accessibility track |
| `class_chart.py` | Bar chart of TEs by class/family |
| `te_family_retroseq.py` | Bar chart of TEs by family (rocket colour scheme) |
| `goterms_socstress_godden.py` | Plotting GO terms analyses from ShinyGO outputs |
| `ppi.py` | Protein-protein interaction plots where TE insertion hits a gene |
| `chrom_rename.py` | Rename chromosomes to numeric Ensembl GRCz11 format |
| `TE_Chrom_count.py` | Count & plot TE insertions per chromosome |
| `features_stats_barchart.py` | Fisher’s exact test on genomic feature overlaps |
| `features_heatmap_2_rows.py` | Heatmap with Fisher’s & Odds ratio testing |
| `retroseq_loci_heatmap.py` | Heatmap of TE insertion loci with enrichment analysis |
| `retroseq_ensembl_vep.py` | Ensembl VEP analysis of Retroseq outputs- heatmap and bar chart |



---

### 🧬 Genome Coverage Analysis  
| Script | Description |
|--------|-------------|
| `average.py`, `coverage.py`, `violinplot.py`, `boxplot.py` | Summarise average depth from `samtools` and plot by chromosome |

---

### 🌱 miRanda Analysis  
| Script | Description |
|--------|-------------|
| `binomial.py` | Calculate enrichment of genes targeted by DE genes & miRNAs |
| `binomialplot.py` | Bar charts of observed vs expected with binomial results |

📖 **Full miRanda analysis** → [See here](https://github.com/alicegodden/paternalsocstress/tree/main/miRanda)  

---

## 📌 Notes  
- Scripts are organised by analysis type for clarity  
- Designed for reproducible zebrafish genomics workflows  
- Some scripts are standalone, others work in pipelines  

---
