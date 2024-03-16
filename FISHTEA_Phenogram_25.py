# Title: Phenogram for Zebrafish chromosomes
# Subtile: Phenogram for overlapping significantly DE genes and TEs
# Author: Dr. Alice M. Godden

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Ellipse
import seaborn as sns

# Set the color palette to Set3
set3_palette = sns.color_palette("tab20")

# Read the matched CSV file
matched_df = pd.read_csv("matched_chromfilter.csv")

# Read the chrom_end.txt file and remove duplicate markers
end_data = []
prev_chrom = None
with open('chrom_end.txt') as f:
    for line in f:
        chrom, pos = line.split()
        if prev_chrom != chrom:
            end_data.append((chrom, int(pos)))
            prev_chrom = chrom

# Create the plot
fig, ax = plt.subplots(figsize=(15, 8))

# Add scatter plot for TE locations
ax.scatter(matched_df['chromosome_TE'], matched_df['TE_start'], label='TE locations', color=set3_palette[0], alpha=1, s=75)

# Add scatter plot for gene locations with a valid color and adjusted transparency
ax.scatter(matched_df['Gene_chromosome'], matched_df['Gene_start'], label='Gene locations', color=set3_palette[1], alpha=0.7, s=25)

# Add centromere markers
centromere_data = pd.read_csv('chrcen.txt', sep='\s+', header=None, names=['chromosome', 'cen_position'])
for _, row in centromere_data.iterrows():
    ax.plot(row['chromosome'], row['cen_position'], marker='D', markersize=5, color=set3_palette[2], zorder=2, alpha=1)

# Set labels and title
plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Position', fontsize=18, fontweight='bold')
plt.title('RNA-seq: Significantly DE overlapping Genes and TEs', fontsize=20, fontweight='bold')

# Set x-axis tick positions and labels dynamically based on unique chromosome values
unique_chromosomes = sorted(matched_df['chromosome_TE'].unique())
ax.set_xticks(range(1, len(unique_chromosomes) + 1))
ax.set_xticklabels(unique_chromosomes, fontsize=14, fontweight='bold')

# Add ellipses for each half of the chromosome
for i, (chrom, end) in enumerate(end_data, start=1):
    centromere_position = centromere_data[centromere_data['chromosome'] == int(chrom)]['cen_position'].iloc[0]

    # Add ellipse for the first half of the chromosome (from 0 to centromere)
    ax.add_patch(
        Ellipse((i, centromere_position / 2), 0.2, centromere_position, edgecolor='grey', linewidth=0.4, fill=False,
                zorder=0))

    # Add ellipse for the second half of the chromosome (from centromere to ylim)
    ax.add_patch(
        Ellipse((i, (end + centromere_position) / 2), 0.2, end - centromere_position, edgecolor='grey', linewidth=0.3,
                fill=False, zorder=0))

# Set y-axis limit to start at 0 and end at the maximum value in the second column of end_data
ax.set_ylim(0)
ax.set_xlim(0.4, len(unique_chromosomes) + 0.6)

# Show the plot
plt.legend()
plt.tight_layout()
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')
plt.savefig("te_genes_dotplot.png", dpi=600)
plt.show()


# Centromere locations chrcen.txt
# 1 33133433
# 2 19743539
# 3 49120158
# 4 25975534
# 5 50391699
# 6 34866988
# 7 60898697
# 8 30099230
# 9 13283057
# 10 9894335
# 11 9507769
# 12 27779269
# 13 19555524
# 14 16933009
# 15 13943685
# 16 19463555
# 17 48225166
# 18 24528805
# 19 19777257
# 20 11319545
# 21 28908825
# 22 21520944
# 23 13631408
# 24 15196971
# 25 20793603

# Chromosome length chrom_end.txt
# 1 59578282
# 2 59640629
# 3 62628489
# 4 78093715
# 5 72500376
# 6 60270059
# 7 74282399
# 8 54304671
# 9 56459846
# 10 45420867
# 11 45484837
# 12 49182954
# 13 52186027
# 14 52660232
# 15 48040578
# 16 55266484
# 17 53461100
# 18 51023478
# 19 48449771
# 20 55201332
# 21 45934066
# 22 39133080
# 23 46223584
# 24 42172926
# 25 37502051

