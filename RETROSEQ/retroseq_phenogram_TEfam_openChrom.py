# Title: Retroseq Phenogram
# Subtitle: With TE family coded points and plotting chromatin accessibility track

# Author: Dr. Alice M. Godden

import matplotlib.pyplot as plt
import pysam
import seaborn as sns
from matplotlib.patches import Ellipse, PathPatch
from matplotlib.lines import Line2D
from matplotlib.path import Path  # Import Path from matplotlib.path
import numpy as np


def extract_chrom_pos_from_vcf(vcf_file):
    chrom_pos_dict = {str(i): [] for i in range(1, 26)}
    te_category_dict = {str(i): [] for i in range(1, 26)}
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            chromosome = record.chrom
            if chromosome in chrom_pos_dict:
                position = record.pos
                info_field = record.info.get('MEINFO', '')
                if isinstance(info_field, tuple):
                    info_field = info_field[0]
                elif not isinstance(info_field, str):
                    info_field = str(info_field)

                if info_field.strip():
                    te_type = info_field.split()[0]
                    if 'DNA' in te_type:
                        te_category = 'DNA'
                    elif 'LINE' in te_type:
                        te_category = 'LINE'
                    elif 'LTR' in te_type:
                        te_category = 'LTR'
                    elif 'SINE' in te_type:
                        te_category = 'SINE'
                    elif 'RC' in te_type:
                        te_category = 'RC'
                    elif 'SATELLITE' in te_type:
                        te_category = 'Satellite'
                    else:
                        te_category = 'Unknown'
                else:
                    te_category = 'Unknown'

                chrom_pos_dict[chromosome].append(position)
                te_category_dict[chromosome].append(te_category)
    return chrom_pos_dict, te_category_dict


def read_dopes_file(dopes_file):
    dopes_data = []
    with open(dopes_file) as f:
        next(f)  # Skip header line
        for line in f:
            chrom, chromStart, chromEnd, name, score, strand = line.split()
            chrom = chrom.replace("chr", "")
            dopes_data.append((chrom, int(chromStart), int(chromEnd), name, int(score), strand))
    return dopes_data


# Replace with your VCF file path
vcf_file = 'MT_HEADERFILTERPY_merged_with_header_USEMEBABY_win150_int_filtered.vcf'
chrom_pos_dict, te_category_dict = extract_chrom_pos_from_vcf(vcf_file)

# Read the chrcen.txt file
chrcen_data = []
with open('chrcen.txt') as f:
    for line in f:
        chrom, pos = line.split()
        chrcen_data.append((chrom, int(pos)))

# Read the chrom_end.txt file
end_data = []
with open('chrom_end.txt') as f:
    for line in f:
        chrom, pos = line.split()
        end_data.append((chrom, int(pos)))

# Color palette for TE types
flare_palette = sns.color_palette("rocket")
te_colors = {
    'DNA': flare_palette[0],
    'LINE': flare_palette[2],
    'LTR': flare_palette[3],
    'RC': flare_palette[4],
    'Satellite': flare_palette[5],
    'SINE': flare_palette[1],
    'Unknown': 'grey'
}

# Create the plot
fig, ax = plt.subplots(figsize=(15, 8))

# Plotting TE insertion points with alpha transparency
for chromosome, positions in chrom_pos_dict.items():
    categories = te_category_dict[chromosome]
    colors = [te_colors[category] for category in categories]
    ax.scatter([int(chromosome)] * len(positions), positions, color=colors, alpha=0.5, label='Chromosome ' + chromosome)

# Add ellipses for each half of the chromosome
ellipse_handles = []
for i, (chrom, end) in enumerate(end_data):
    centromere_position = next(cen for chr_name, cen in chrcen_data if chr_name == chrom)

    # Add ellipse for the first half of the chromosome (from 0 to centromere)
    ax.add_patch(
        Ellipse((int(chrom), centromere_position / 2), 0.27, centromere_position, edgecolor='black', linewidth=0.4,
                fill=False,
                zorder=20))

    # Add ellipse for the second half of the chromosome (from centromere to end)
    ax.add_patch(
        Ellipse((int(chrom), (end + centromere_position) / 2), 0.27, end - centromere_position, edgecolor='black',
                linewidth=0.4,
                fill=False, zorder=20))

# Add centromere diamonds
for chrom, cen in chrcen_data:
    ax.plot([int(chrom)], [cen], marker='D', markersize=3.75, markerfacecolor="grey", color='grey', zorder=21)

# Add chromatin accessibility data next to the ellipse
for chrom, chromStart, chromEnd, name, score, strand in read_dopes_file('daniocode_hub_280355_dopes_all.txt'):
    ax.hlines(y=[chromStart, chromEnd], xmin=int(chrom) + 0.15, xmax=int(chrom) + 0.25, color='darkred', linewidth=0.2)

# Customizing plot details
plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Position', fontsize=18, fontweight='bold')
plt.title('Overlay Plot: Non-reference TE Insertion Locations vs. Chromatin Accessibility', fontsize=20,
          fontweight='bold')
plt.xticks(range(1, 26), [chrom for chrom, _ in end_data], fontsize=18, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')

# Customize legend for both plots
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, alpha=0.5, label=te_type)
    for te_type, color in te_colors.items() if te_type != 'Unknown'
] + [
    Line2D([0], [0], color='darkred', lw=2, label='Open regions'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='grey', markersize=8, label='Centromere')
]
plt.legend(handles=legend_elements, fontsize=12, loc='upper right')

# Save and display the plot
output_file = 'overlay_plot_with_centromere_alpha.png'  # Replace with your desired filename and extension
plt.savefig(output_file, dpi=600)
plt.show()

# Note: chromatin accessibility track:
# chromatin track was accessed from DANIO code hub_280355_dopes_all on UCSC browser
# Unmarked open chromatin regions
# This collection contains two tracks with regions which have an ATAC-seq signal, but without observable CRE-associated chromatin marks.
#
# COPEs: Constitutive Orphan Predicted Elements (COPEs) are regions which are constitutively open throughout development.
#
# pooled DOPEs: Dynamic Orphan Predicted Elements (DOPEs) are regions which are open only in specific developmental stages. DOPEs from all stages were pooled together into this track.
#
# For more information about the origin of the data, see our publication: Baranasic, Damir, et al. "Integrated annotation and analysis of genomic features reveal new types of functional elements and large-scale epigenetic phenomena in the developing zebrafish."
