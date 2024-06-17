import matplotlib.pyplot as plt
import pysam
import seaborn as sns
from matplotlib.patches import Ellipse

def extract_chrom_pos_from_vcf(vcf_file):
    chrom_pos_dict = {str(i): [] for i in range(1, 26)}
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            chromosome = record.chrom
            if chromosome in chrom_pos_dict:
                position = record.pos
                chrom_pos_dict[chromosome].append(position)
    return chrom_pos_dict

def get_chrom_end(chrom):
    with open('chrom_end.txt') as f:
        for line in f:
            chrom_name, end_pos = line.split()
            if chrom_name == chrom:
                return int(end_pos)
    return None

# Replace 'your_vcf_file.vcf.gz' with the path to your VCF file or .vcf.gz
vcf_file = 'high_filtered_sift_CC_snpeff_female_allchrs.eff.vcf'
chrom_pos_dict = extract_chrom_pos_from_vcf(vcf_file)

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

# Color palette for chromosomes
palette = sns.color_palette("rocket", n_colors=len(chrom_pos_dict))

# Create the plot
fig, ax = plt.subplots(figsize=(15, 8))

# Plot chromosome positions with unique colors
for i, (chromosome, positions) in enumerate(chrom_pos_dict.items()):
    color = palette[i]
    plt.scatter([int(chromosome)] * len(positions), positions, color=color, label='Chromosome ' + chromosome)

# Add ellipses for chromosome halves
for i, (chrom, end) in enumerate(end_data):
    centromere_position = next(cen for chr_name, cen in chrcen_data if chr_name == chrom)

    # Add ellipse for the first half of the chromosome (from 0 to centromere)
    ax.add_patch(
        Ellipse((int(chrom), centromere_position / 2), 0.27, centromere_position, edgecolor='black', linewidth=0.4,
                fill=False, zorder=20)
    )

    # Add ellipse for the second half of the chromosome (from centromere to end)
    ax.add_patch(
        Ellipse((int(chrom), (end + centromere_position) / 2), 0.27, end - centromere_position, edgecolor='black',
                linewidth=0.4, fill=False, zorder=20)
    )

# Add centromere diamonds
for chrom, cen in chrcen_data:
    ax.plot([int(chrom)], [cen], marker='D', markersize=5, markerfacecolor="grey", color='grey', zorder=21)

# Add genes data as text labels with marks
for chrom, pos, gene in genes_data:
    x, y = int(chrom), pos
    ax.annotate(gene, xy=(x, y), ha='center', va='center', fontsize=11, fontweight='bold', zorder=22,
                arrowprops=dict(arrowstyle="<-", color='red', facecolor='none', linewidth=0))

# Customize plot details
plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Position', fontsize=18, fontweight='bold')
plt.title('SNPEFF & SNPSIFT High Impact = < 0.05 - Ovaries temperature vs control', fontsize=20, fontweight='bold')
plt.plot([], [], label='', color='none')  # Remove legend by setting an empty label
plt.tight_layout()
plt.xticks(range(1, 26), [chrom for chrom, _ in end_data], fontsize=18, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')

# Save and display the plot
output_file = 'female_highimpact_snpeffsift_with_ellipses_and_centromeresrocket.png'
plt.savefig(output_file)
plt.show()
