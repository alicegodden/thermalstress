# Bash scripts
TBC


# Environmental stress scripts
Thermal stress Zebrafish genomics project


# Social stress in Zebrafish transcriptome and sperm sRNAs - R notebook
Godden_seq.Rmd
ADD PHENOGRAM R SCRIPT

# Thermal stress in Zebrafish gonads - R notebook for RNA-seq and sRNA analyses
AG_RNAseq.Rmd


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



# RETROSEQ analysis

filter_retroseq.py # filtering on GQ and FL Retroseq outputs to filter vcf for most confident novel TE insertion calls

retroseq.py # Plotting Retroseq outputs

class_chart.py # Plotting bar chart of TEs by class

goterms_socstress_godden.py # Plotting GO terms analyses from ShinyGo outputs

ppi.py # plotting Retroseq insertion protein-protein interaction terms, where mutation hits a gene


GENOME COVERAGE analysis

average.py , coverage.py, violinplot.py and boxplot.py # For summarising the average depth from samtools and making plots to analyse this by chromosome at whole genome level
