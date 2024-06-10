# Bash scripts
TBC


# Environmental stress scripts
Thermal stress Zebrafish genomics project


# Social stress in Zebrafish transcriptome and sperm sRNAs - R notebook
Godden_seq.Rmd


# Thermal stress in Zebrafish gonads - R notebook for RNA-seq and sRNA analyses
AG_RNAseq.Rmd
telescopeCounts.R - Telescope output readin to make universal counts file across multiple samples


# Thermal stress in Zebrafish gonads - R notebook- for DNA-seq analyses of variants
Haploid analyses.Rmd


# Python scripts 

# SNPEFF analysis

filter.py # Filtering SnpEff and SnpSift outputs

tstvplot.py # Plotting Ts/Tv ratios output from SnpEff

snpeff.py # Plotting variants from SnpEff and SnpSift outputs

autobubble_goplot.py # Plotting go terms outputs (.csv from ShinyGO) with annotated p-values


# FishTEA- Fish Transposable element analyzer scripts

FISHTEA_annotate_family.py # Annotating the significantly differentially expressed TE families from TEtranscripts pipeline with chromosome and genomic co-ordinates

FISHTEA_matching.py # Matching the TEs with genes that are overlapping in location

FISHTEA_Phenogram_25.py # Plotting overlapping significantly DE genes and TEs on Zebrafish genome

FISHTEA_chromatin_phenogram.py # plotting open/accessible regions of the genome 



# RETROSEQ analysis

filter_retroseq.py # filtering on GQ and FL Retroseq outputs to filter vcf for most confident novel TE insertion calls

retroseq.py and vcfinput_phenogram.py (for chromosome shapes) # Plotting Retroseq outputs

class_chart.py # Plotting bar chart of TEs by class

goterms_socstress_godden.py # Plotting GO terms analyses from ShinyGo outputs

ppi.py # plotting Retroseq insertion protein-protein interaction terms, where mutation hits a gene

chrom_rename.py # For renaming chromosomes to numeric IDs to match Ensembl formatting GRCz11 genome


GENOME COVERAGE analysis

average.py , coverage.py, violinplot.py and boxplot.py # For summarising the average depth from samtools and making plots to analyse this by chromosome at whole genome level

# Polar bear analysis

bearphenogram.py # Plotting a basic phenogram of a list of significantly differentially expressed TEs from Telescope outputs, genome ASM1731132v1

chrom_end_bear.txt # Chromosome/scaffold legnths text file, genome ASM1731132v1

bear_te_class.py # Plotting bar charts of TEs


