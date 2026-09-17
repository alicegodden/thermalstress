import matplotlib.pyplot as plt
import pysam
import seaborn as sns

def extract_chrom_pos_from_vcf(vcf_file):
    chrom_pos_dict = {str(i): [] for i in range(1, 26)}
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            chromosome = record.chrom
            if chromosome in chrom_pos_dict:
                position = record.pos
                chrom_pos_dict[chromosome].append(position)
    return chrom_pos_dict

def read_chrom_lengths(chrom_file):
    """Read chromosome lengths from file and return a dictionary with lengths in Mb."""
    chrom_lengths = {}
    with open(chrom_file, 'r') as f:
        for line in f:
            chrom, length = line.strip().split()
            chrom_lengths[chrom] = int(length) / 1_000_000  # Convert to Mb
    return chrom_lengths

# Replace with your VCF file paths
vcf_files = [
    'ctrl_female_exp_unique_8_win100_gq750_fl8.vcf',  # Replace with your file paths
    'female_exp_unique_8_win100_filterpy_gq750_fl8.vcf',
    'ctrl_male_exp_unique_8_win100_gq750_fl8.vcf',
    'male_exp_unique_8_win100_filterpy_gq750_fl8.vcf'
]
# Corresponding sample names
sample_names = [
    'Control_Ovaries',
    'Thermal_Ovaries',
    'Control_Testes',
    'Thermal_Testes'
]

# Load chromosome lengths (assuming chrom_end.txt has chromosome lengths)
chrom_lengths = read_chrom_lengths('chrom_end.txt')

# Store normalized TE counts for each VCF file
all_te_counts = []

# Process each VCF file
for vcf_file in vcf_files:
    chrom_pos_dict = extract_chrom_pos_from_vcf(vcf_file)
    te_counts = {chrom: len(positions) / chrom_lengths.get(chrom, 1) for chrom, positions in chrom_pos_dict.items()}  # Normalize by chromosome length
    all_te_counts.append(te_counts)

# Prepare data for plotting
chromosomes = list(map(str, range(1, 26)))
data = {sample_name: [te_counts.get(chrom, 0) for chrom in chromosomes] for sample_name, te_counts in zip(sample_names, all_te_counts)}

# Color palette for the lines
rocket_palette = sns.color_palette("rocket", len(sample_names))

# Plotting
fig, ax = plt.subplots(figsize=(15, 8))

for sample_name, color in zip(sample_names, rocket_palette):
    ax.plot(chromosomes, data[sample_name], label=sample_name, color=color, marker='o')

# Customizing plot details
plt.xlabel('Chromosome', fontsize=18, fontweight='bold')
plt.ylabel('TE counts per Mb', fontsize=18, fontweight='bold')  # Updated label
plt.title('TEs per Chromosome (Normalized by Length) from Retroseq', fontsize=20, fontweight='bold')
plt.xticks(fontsize=12, fontweight='bold')
plt.yticks(fontsize=12, fontweight='bold')
plt.legend(title='Samples', fontsize=12, title_fontsize=14, loc='upper right')

# Save and display the plot
output_file = 'TE_counts_per_chromosome_per_MB_fish_HS_retroseq.png'
plt.savefig(output_file, dpi=1200)
plt.show()
