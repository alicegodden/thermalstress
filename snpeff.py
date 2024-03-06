# Title- snpeff.py
# Author - Dr. Alice Godden
import matplotlib.pyplot as plt
import pysam

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
    x_end, y_end = chrom, pos
    # ax.plot([x, x_end], [y, y_end], color='blue')
    ax.annotate(gene, xy=(x, y), ha='center', va='center', fontsize=11.5, fontweight='bold',
                arrowprops=dict(arrowstyle="<-", color='red', facecolor='none', linewidth=0))

# Add the chromend data as marks on the plot
for chrom, end in end_data:
    end_pos = get_chrom_end(chrom)
    if end_pos is not None:
        plt.plot([chrom], [end_pos], marker='_', markersize=10, color='black')

plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('Position', fontsize=18, fontweight='bold')
plt.title('SNPEFF & SNPSIFT High confidence p= < 0.05 - Ovaries temperature vs control', fontsize=20, fontweight='bold')
# Remove the legend by setting an empty label
plt.plot([], [], label='', color='none')  # Adjust the color as needed
# #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5)) # to have a legend of the chromosome colors on the right
plt.tight_layout()  # Adjust the layout to avoid overlapping
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')


# save the file
output_file = 'female_high_filt_sig_snpsift_snpeff.png'  # Replace with your desired filename and extension
plt.savefig(output_file)
plt.show()
