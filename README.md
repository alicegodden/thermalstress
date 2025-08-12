
# Environmental stress scripts
Thermal stress Zebrafish genomics project

# Contents

# Water analysis
water_testing.py # water testing parameters and graphs


# /Rscripts
AG_RNAseq.Rmd # R notebook for DESeq2 analyses

manhattan_99.R # Manhattan plots with 99th percentile plotting

pheatmap_logcpm_z_heatstress.R # Normalised counts and heatmaps from DESeq2 outputs

tajD_fst_genes.R # Tajima's D and FST plots tracks with genes at chromosomal level

telescopeCounts.R # Custome scripts for analysis and assimilation of count data from Telescope pipeline

testing_glm_HS_retroseq.R # GLM of Retroseq data at family level


# Thermal stress in Zebrafish gonads - R notebook- for DNA-seq analyses of variants
Haploid analyses.Rmd


# Python scripts 

# SNPEFF analysis

filter.py # Filtering SnpEff and SnpSift outputs

tstvplot.py # Plotting Ts/Tv ratios output from SnpEff

snpeff.py # Plotting variants from SnpEff and SnpSift outputs

vcfinput_phenogram.py # Phenogram plot

autobubble_goplot.py # Plotting go terms outputs (.csv from ShinyGO) with annotated p-values

unique.py # filtering non-unique rows from experimental retroseq outputs to generate file of unique variants


# FishTEA- Fish Transposable element analyzer scripts

FISHTEA_annotate_family.py # Annotating the significantly differentially expressed TE families from TEtranscripts pipeline with chromosome and genomic co-ordinates

FISHTEA_matching.py # Matching the TEs with genes that are overlapping in location

FISHTEA_Phenogram_25.py # Plotting overlapping significantly DE genes and TEs on Zebrafish genome

FISHTEA_chromatin_phenogram.py # plotting open/accessible regions of the genome 



# RETROSEQ analysis

filter_retroseq.py # filtering on GQ and FL Retroseq outputs to filter vcf for most confident novel TE insertion calls

retroseq.py and vcfinput_phenogram.py (for chromosome shapes) # Plotting Retroseq outputs

retroseq_phenogram_TEcolours.py # Plotting Retroseq outputs and grouping TE's by family, plotting points in colours based on family

retroseq_phenogram_TEfam_openChrom.py # As above but with the addition of a chromatin accessibility track

class_chart.py # Plotting bar chart of TEs by class/family

te_family_retroseq.py # Plotting bar chart of TEs by family, rocket colour scheme

goterms_socstress_godden.py # Plotting GO terms analyses from ShinyGo outputs

ppi.py # plotting Retroseq insertion protein-protein interaction terms, where mutation hits a gene

chrom_rename.py # For renaming chromosomes to numeric IDs to match Ensembl formatting GRCz11 genome

TE_Chrom_count.py # Counting and plotting number of non-reference TE insertions mutations per chromosome

features_stats_barchart.py # Fisher's exact testing on overlapping genomic features of TE non-reference insertion mutations

features_heatmap_2_rows.py # Heatmap with Fishers and Odds ratio resting for enrichment of features that TE insertion mutations overlap


GENOME COVERAGE analysis

average.py , coverage.py, violinplot.py and boxplot.py # For summarising the average depth from samtools and making plots to analyse this by chromosome at whole genome level

# miRanda analysis
binomial.py # calculate enrichment of genes targeted by our significantly differentially expressed genes and miRNAs

binomialplot.py # plotting bar charts observed v expected and binomial results

Full miRanda analysis as described previously [here](https://github.com/alicegodden/paternalsocstress/tree/main/miRanda)




