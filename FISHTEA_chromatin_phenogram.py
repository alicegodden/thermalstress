# Title: Chromatin phenogram
# Subtitle: DanRer11 genome accessibility throughout development
# Author: Dr. Alice M. Godden

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D


def read_dopes_file(dopes_file):
    dopes_data = []
    with open(dopes_file) as f:
        next(f)  # Skip header line
        for line in f:
            chrom, chromStart, chromEnd, name, score, strand = line.split()
            # Remove the "chr" prefix to match the chromosome names in other files
            chrom = chrom.replace("chr", "")
            dopes_data.append((chrom, int(chromStart), int(chromEnd), name, int(score), strand))
    return dopes_data

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

# Read the dopes.txt file
dopes_file = 'daniocode_hub_280355_dopes_all.txt'
dopes_data = read_dopes_file(dopes_file)

# Create the plot
fig, ax = plt.subplots(figsize=(15, 8))

# Add the dopes data as horizontal lines on the plot
for chrom, chromStart, chromEnd, name, score, strand in dopes_data:
    ax.hlines(y=[chromStart, chromEnd], xmin=int(chrom) - 0.08, xmax=int(chrom) + 0.08, color='darkred', linewidth=0.5)

# Add the chrcen data as marks on the plot
for chrom, cen in chrcen_data:
    ax.plot([int(chrom)], [cen], marker='D', markersize=5, markerfacecolor="grey", color='black', zorder=21)

# Add ellipses for each half of the chromosome
for i, (chrom, end) in enumerate(end_data):
    centromere_position = next(cen for chr_name, cen in chrcen_data if chr_name == chrom)

    # Add ellipse for the first half of the chromosome (from 0 to centromere)
    ax.add_patch(
        Ellipse((int(chrom), centromere_position / 2), 0.27, centromere_position, edgecolor='black', linewidth=0.4, fill=False,
                zorder=20))

    # Add ellipse for the second half of the chromosome (from centromere to end)
    ax.add_patch(
        Ellipse((int(chrom), (end + centromere_position) / 2), 0.27, end - centromere_position, edgecolor='black', linewidth=0.4,
                fill=False, zorder=20))

plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Position', fontsize=18, fontweight='bold')
plt.title('Open regions of chromatin danRer11', fontsize=20, fontweight='bold')
plt.tight_layout()
plt.xticks(range(1, 26), [chrom for chrom, _ in end_data], fontsize=18, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')

# Customize legend
legend_elements = [
    Line2D([0], [0], color='darkred', lw=2, label='Open regions'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='grey', markersize=8, label='Centromere')
]
plt.legend(handles=legend_elements, fontsize=16, loc='upper right')


# Save the file
output_file = 'chromatin_phenogram_danio.png'  # Replace with your desired filename and extension
plt.savefig(output_file)
plt.show()


# chromatin track was accessed from DANIO code hub_280355_dopes_all on UCSC browser
# Unmarked open chromatin regions
# This collection contains two tracks with regions which have an ATAC-seq signal, but without observable CRE-associated chromatin marks.
#
# COPEs: Constitutive Orphan Predicted Elements (COPEs) are regions which are constitutively open throughout development.
#
# pooled DOPEs: Dynamic Orphan Predicted Elements (DOPEs) are regions which are open only in specific developmental stages. DOPEs from all stages were pooled together into this track.
#
# For more information about the origin of the data, see our publication: Baranasic, Damir, et al. "Integrated annotation and analysis of genomic features reveal new types of functional elements and large-scale epigenetic phenomena in the developing zebrafish."
