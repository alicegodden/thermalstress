# Title: Retroseq phenogram
# Subtitle : TE family colour coding

# Author: Dr. Alice M. Godden

import matplotlib.pyplot as plt
import pysam
import seaborn as sns
from matplotlib.patches import Ellipse


def extract_chrom_pos_from_vcf(vcf_file):
    chrom_pos_dict = {str(i): [] for i in range(1, 26)}
    te_category_dict = {str(i): [] for i in range(1, 26)}
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            chromosome = record.chrom
            if chromosome in chrom_pos_dict:
                position = record.pos
                info_field = record.info.get('MEINFO', '')  # Ensure info_field is handled as expected

                # Handle different types of info_field
                if isinstance(info_field, tuple):
                    info_field = info_field[0]  # Assuming it's a tuple, take the first element
                elif not isinstance(info_field, str):
                    info_field = str(info_field)  # Convert to string if not already

                # Extract only the first term from MEINFO
                if info_field.strip():  # Check if info_field is not empty
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


# Replace with your VCF file path
vcf_file = 'FT_HEADERFILTERPY_merged_with_header_USEMEBABY_win150_int_filtered.vcf'
chrom_pos_dict, te_category_dict = extract_chrom_pos_from_vcf(vcf_file)

# Read the chrcen.txt file
chrcen_data = []
with open('chrcen.txt') as f:
    for line in f:
        chrom, pos = line.split()
        chrcen_data.append((chrom, int(pos)))

# Read the genes.txt file
genes_data = []
with open('genes.txt') as f:
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
    'SINE': flare_palette[1]
    # Remove 'Unknown' from here
}

# Create the plot
fig, ax = plt.subplots(figsize=(15, 8))
for chromosome, positions in chrom_pos_dict.items():
    categories = te_category_dict[chromosome]
    colors = [te_colors[category] for category in categories]
    plt.scatter([chromosome] * len(positions), positions, color=colors, label='Chromosome' + chromosome)

# Add the chrcen data as marks on the plot
for chrom, cen in chrcen_data:
    ax.plot([chrom], [cen], marker='D', markersize=5, markerfacecolor="grey", color='black')

# Add the genes data as text labels with marks
for chrom, pos, gene in genes_data:
    x, y = chrom, pos
    ax.annotate(gene, xy=(x, y), ha='center', va='center', fontsize=11.5, fontweight='bold',
                arrowprops=dict(arrowstyle="<-", color='red', facecolor='none', linewidth=0))

# Add ellipses for each half of the chromosome
for i, (chrom, end) in enumerate(end_data, start=0):
    centromere_position = chrcen_data[i][1]  # Assuming centromere_data aligns with end_data

    # Add ellipse for the first half of the chromosome (from 0 to centromere)
    ax.add_patch(
        Ellipse((i, centromere_position / 2), 0.2, centromere_position, edgecolor='grey', linewidth=0.4, fill=False,
                zorder=0))

    # Add ellipse for the second half of the chromosome (from centromere to ylim)
    ax.add_patch(
        Ellipse((i, (end + centromere_position) / 2), 0.2, end - centromere_position, edgecolor='grey', linewidth=0.3,
                fill=False, zorder=0))

plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Position', fontsize=18, fontweight='bold')
plt.title('Non-reference TE insertion Locations Ovaries Temperature', fontsize=20, fontweight='bold')
plt.tight_layout()
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')

# Add legend (excluding 'Unknown')
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=te_type)
           for te_type, color in te_colors.items() if te_type != 'Unknown']
ax.legend(handles=handles, title="TE Family", fontsize=12, title_fontsize=14, loc='upper right')

# Save the file
output_file = 'Retroseq_Ovaries_Temperature_phenogram_insertions_colors.png'  # Replace with your desired filename and extension
plt.savefig(output_file)
plt.show()
