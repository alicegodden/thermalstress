# Title : Retroseq chromosomal plot with hypergeometric testing
# Author : Dr. Alice M. Godden

import matplotlib.pyplot as plt
import pysam
import seaborn as sns
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
from scipy.stats import hypergeom
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

def perform_hypergeometric_test(observed, total_genes, expected, total_population):
    """Perform a hypergeometric test for the observed and expected counts."""
    p_value = hypergeom.sf(observed - 1, total_population, expected, total_genes)
    return p_value

# Replace with your VCF file paths
exp_vcf_file = 'female_exp_unique_8_win100_filterpy_gq750_fl8.vcf'
ctrl_vcf_file = 'ctrl_female_exp_unique_8_win100_gq750_fl8.vcf'  # Replace with the actual control VCF file path

# Extract chromosome positions from VCF
exp_chrom_pos_dict, exp_te_category_dict = extract_chrom_pos_from_vcf(exp_vcf_file)
ctrl_chrom_pos_dict, _ = extract_chrom_pos_from_vcf(ctrl_vcf_file)

# Read the chrcen.txt file
chrcen_data = []
with open('chrcen.txt') as f:
    for line in f:
        chrom, pos = line.split()
        chrcen_data.append((chrom, int(pos)))

 #Read the genes.txt file
genes_data = []
with open('genes_poster') as f:
     for line in f:
         chrom, pos, gene = line.split()
         genes_data.append((chrom, int(pos), gene))

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

# Total control count (based on the control VCF)
total_population = sum(len(positions) for positions in ctrl_chrom_pos_dict.values())
total_successes = sum(len(positions) for positions in exp_chrom_pos_dict.values())  # Total TE insertions in exp

# Calculate enriched regions for experimental insertions compared to control
enriched_regions = []
window_size = 25000  # window size 40kb
step_size = 50000  # step size 10kb

# Calculate enriched regions for each chromosome
for chromosome, positions in exp_chrom_pos_dict.items():
    if chromosome not in ctrl_chrom_pos_dict:
        continue  # Skip chromosomes that are not in control

    ctrl_positions = ctrl_chrom_pos_dict[chromosome]
    max_position = end_data[int(chromosome) - 1][1]  # Get the max length of the chromosome directly
    positions.sort()

    # Slide window across the chromosome
    for start in range(0, max_position, step_size):
        end = min(start + window_size, max_position)
        count_input_genes = len([pos for pos in positions if start <= pos < end])
        count_all_genes = len([pos for pos in ctrl_positions if start <= pos < end])

        # Perform hypergeometric test
        p_value = perform_hypergeometric_test(count_input_genes, len(positions), count_all_genes, total_population)

        # If the region is statistically significant (p-value < 0.05), add it to the enriched regions list
        if p_value < 0.01:
            enriched_regions.append((chromosome, start, end, p_value))

# Plotting TE insertion points for experimental data with colors
for chromosome, positions in exp_chrom_pos_dict.items():
    categories = exp_te_category_dict[chromosome]
    colors = [te_colors[category] if category in te_colors else te_colors['Unknown'] for category in categories]
    ax.scatter([int(chromosome)] * len(positions), positions, color=colors, alpha=0.5)

# Highlight significant enriched regions next to the chromosome
for chrom, start, end, p_value in enriched_regions:
    ax.bar(int(chrom) - 0.2, end - start, bottom=start, color='grey', alpha=1, width=0.4,
           label='Enriched Region' if chrom == enriched_regions[0][0] and start == enriched_regions[0][1] else "")

# Add ellipses for each half of the chromosome with space below
for i, (chrom, end) in enumerate(end_data):
    centromere_position = next(cen for chr_name, cen in chrcen_data if chr_name == chrom)

    # Add ellipse for the first half of the chromosome
    ax.add_patch(
        Ellipse((int(chrom), centromere_position / 2), 0.27, centromere_position, edgecolor='black', linewidth=0.6,
                fill=False, zorder=20))

    # Add ellipse for the second half of the chromosome
    ax.add_patch(
        Ellipse((int(chrom), (end + centromere_position) / 2), 0.27, end - centromere_position, edgecolor='black',
                linewidth=0.6, fill=False, zorder=20))


# Add centromere diamonds
for chrom, cen in chrcen_data:
    ax.plot([int(chrom)], [cen], marker='D', markersize=3.75, markerfacecolor="grey", color='grey', zorder=21)

# Add chromatin accessibility data
for chrom, chromStart, chromEnd, name, score, strand in read_dopes_file('daniocode_hub_280355_dopes_all.txt'):
    ax.hlines(y=[chromStart, chromEnd], xmin=int(chrom) + 0.15, xmax=int(chrom) + 0.25, color='darkred', linewidth=0.2)

# After plotting your data and before displaying the plot

# Add ellipses for each half of the chromosome with space below
for i, (chrom, end) in enumerate(end_data):
    centromere_position = next(cen for chr_name, cen in chrcen_data if chr_name == chrom)

    # Add ellipse for the first half of the chromosome
    ax.add_patch(
        Ellipse((int(chrom), centromere_position / 2), 0.27, centromere_position, edgecolor='black', linewidth=0.6,
                fill=False, zorder=20))

    # Add ellipse for the second half of the chromosome
    ax.add_patch(
        Ellipse((int(chrom), (end + centromere_position) / 2), 0.27, end - centromere_position, edgecolor='black',
                linewidth=0.6, fill=False, zorder=20))

# Add the genes data as text labels with marks
for chrom, pos, gene in genes_data:
    x, y = int(chrom) - 0.22, pos  # Shift x position to the left of the ellipse
    ax.annotate(
        gene,
        xy=(int(chrom), pos),  # Use the original gene position for the arrow
        ha='center',
        va='center',
        fontsize=8,
        fontweight='bold',
        fontstyle='italic',
        rotation=90,  # Rotate the label 90 degrees
        xytext=(-15, 0),  # Move the text to the left of the gene position
        textcoords='offset points',  # Use offset points for the text position
        arrowprops=dict(arrowstyle="simple", color='black', linewidth=0.5),  # Ensure arrow is visible
    )

# Adjust y-axis limits to provide additional space below
y_min = min(min(positions) for positions in exp_chrom_pos_dict.values()) - 4000000  # Set to a negative value to create space below ellipses
ax.set_ylim(y_min)  # Set the y-axis limits

# Customizing plot details
plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Position', fontsize=18, fontweight='bold')
plt.title('Non-reference TE Insertion mutations with Enriched Regions', fontsize=20, fontweight='bold')
plt.xticks(range(1, 26), [chrom for chrom, _ in end_data], fontsize=18, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')

# Customize legend for both plots
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, alpha=0.5, label=te_type)
    for te_type, color in te_colors.items() if te_type != 'Unknown'
] + [
    Line2D([0], [0], color='grey', lw=2, label='Enriched Regions'),
    Line2D([0], [0], color='darkred', lw=2, label='Open Chromatin'),
    Line2D([0], [0], marker='D', color='w', markerfacecolor='grey', markersize=8, label='Centromere')
]
plt.legend(handles=legend_elements, fontsize=12, loc='upper right')

# Save and display the plot
output_file = 'test.png'  # Replace with your desired filename and extension
plt.savefig(output_file, dpi=600)

plt.show()
