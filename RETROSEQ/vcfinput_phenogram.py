# Title: Phenogram plotting vcf file inputs
# Author: Dr. Alice M. Godden
# Import libraries
import matplotlib.pyplot as plt
import pysam
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

# Replace 'your_vcf_file.vcf' with the path to your VCF file
vcf_file = 'high_filtered_sift_CC_snpeff_male_allchrs.eff.vcf'
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

# Create the plot
fig, ax = plt.subplots(figsize=(15, 8))
for chromosome, positions in chrom_pos_dict.items():
    plt.scatter([chromosome] * len(positions), positions, label='Chromosome' + chromosome)

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
plt.title('SNPEFF & SNPSIFT High Confidence p = < 0.05 - Testes temperature vs control', fontsize=20, fontweight='bold')
plt.tight_layout()
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')

# Save the file
output_file = 'snpsift_testes_TvC_sig.png'  # Replace with your desired filename and extension
plt.savefig(output_file)
plt.show()
