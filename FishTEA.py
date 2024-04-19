# Title: FishTEA
# Subtitle: Fish transposable element analyzer
# Author: Dr. Alice M. Godden

# Step 1 of 5
import pandas as pd

# Define file paths
sigTEs_file = "sigTEs_positions_telescope_family_TESTES.csv"
sigGenes_file = "te_genes.csv"
log_file = "matching_log_ov.txt"

# Read the CSV files
print("Reading CSV files...")
sigTEs_df = pd.read_csv(sigTEs_file)
sigGenes_df = pd.read_csv(sigGenes_file)

# Create an empty list to store the matched TEs and genes
matched_records = []

# Iterate over each row in sigTEs_df
print("Matching TEs with genes...")
for index_TE, TE_row in sigTEs_df.iterrows():
    TE_start = TE_row['TE_start']
    TE_end = TE_row['TE_end']
    TE_chromosome = TE_row['chromosome_TE']

    # Iterate over each row in sigGenes_df
    for index_gene, gene_row in sigGenes_df.iterrows():
        gene_start = gene_row['Gene_start']
        gene_end = gene_row['Gene_end']
        gene_chromosome = gene_row['Gene_chromosome']

        # Check if the chromosome information matches
        if TE_chromosome == gene_chromosome:
            # Check for overlap between TE and gene
            if (TE_start <= gene_end) and (gene_start <= TE_end):
                # Store the matched TE and gene information
                matched_records.append({
                    'chromosome_TE': TE_chromosome,
                    'TE_name': TE_row['TE_name'],
                    'TE_class': TE_row['TE_class'],
                    'TE_family': TE_row['TE_family'],
                    'TE_start': TE_start,
                    'TE_end': TE_end,
                    'Gene_ID': gene_row['Gene_ID'],
                    'Gene_name': gene_row['Gene_name'],
                    'Gene_start': gene_start,
                    'Gene_end': gene_end,
                    'Gene_chromosome': gene_chromosome
                })

    # Print progress update
    print(f"Processed {index_TE + 1} out of {len(sigTEs_df)} TEs.", end="\r")

# Create a DataFrame from the matched_records list
matched_df = pd.DataFrame(matched_records)

# Check if any matches were found
if len(matched_df) == 0:
    print("No TEs and genes are overlapping in position, sorry")
else:
    # Write the matched dataframe to a new CSV file
    matched_df.to_csv("matched_TESTES.csv", index=False)

    # Write log information to a log file
    with open(log_file, "w") as f:
        f.write("Matching TEs with genes completed.\n")
        f.write(f"Number of TEs processed: {len(sigTEs_df)}\n")
        f.write(f"Number of matched TEs and genes: {len(matched_df)}\n")

    print("\nMatching completed. Results saved to 'matched_TESTES.csv'.")


# Step 2 of 5
# Matching overlaping sig DE genes and sig DE TEs

import pandas as pd

# Define file paths
sigTEs_file = "sigTEs_positions.csv"
sigGenes_file = "sigGenes.csv"
log_file = "matching_log.txt"

# Read the CSV files
print("Reading CSV files...")
sigTEs_df = pd.read_csv(sigTEs_file)
sigGenes_df = pd.read_csv(sigGenes_file)

# Create an empty list to store the matched TEs and genes
matched_records = []

# Iterate over each row in sigTEs_df
print("Matching TEs with genes...")
for index_TE, TE_row in sigTEs_df.iterrows():
    TE_start = TE_row['TE_start']
    TE_end = TE_row['TE_end']
    TE_chromosome = TE_row['chromosome_TE']

    # Iterate over each row in sigGenes_df
    for index_gene, gene_row in sigGenes_df.iterrows():
        gene_start = gene_row['Gene_start']
        gene_end = gene_row['Gene_end']
        gene_chromosome = gene_row['Gene_chromosome']

        # Check if the chromosome information matches
        if TE_chromosome == gene_chromosome:
            # Check for overlap between TE and gene
            if (TE_start <= gene_end) and (gene_start <= TE_end):
                # Store the matched TE and gene information
                matched_records.append({
                    'chromosome_TE': TE_chromosome,
                    'TE_name': TE_row['TE_name'],
                    'TE_class': TE_row['TE_class'],
                    'TE_family': TE_row['TE_family'],
                    'TE_start': TE_start,
                    'TE_end': TE_end,
                    'Gene_ID': gene_row['Gene_ID'],
                    'Gene_name': gene_row['Gene_name'],
                    'Gene_start': gene_start,
                    'Gene_end': gene_end,
                    'Gene_chromosome': gene_chromosome
                })

    # Print progress update
    print(f"Processed {index_TE + 1} out of {len(sigTEs_df)} TEs.", end="\r")

# Create a DataFrame from the matched_records list
matched_df = pd.DataFrame(matched_records)

# Write the matched dataframe to a new CSV file
matched_df.to_csv("matched.csv", index=False)

# Write log information to a log file
with open(log_file, "w") as f:
    f.write("Matching TEs with genes completed.\n")
    f.write(f"Number of TEs processed: {len(sigTEs_df)}\n")
    f.write(f"Number of matched TEs and genes: {len(matched_df)}\n")

print("\nMatching completed. Results saved to 'matched_TEchrom.csv'.")

# Step 3 of 5
# Phenogram plotting results
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Ellipse
import seaborn as sns

# Set the colorblind-friendly color palette
colorblind_palette = sns.color_palette("colorblind")
diverging_colors = sns.color_palette("RdBu", 10)
sns.palplot(diverging_colors)

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

# Set colors based on TE_family column
colors = []
for family in matched_df['TE_family']:
    if 'SINE' in family or 'LINE' in family or 'LTR' in family:
        colors.append(diverging_colors[0])  # Blue
    elif 'DNA' in family or 'RC' in family:
        colors.append(diverging_colors[6])  # Red
    else:
        colors.append(diverging_colors[3])  # Green (default)

# Add scatter plot for TE locations with different colors
# Create a list to store marker sizes based on TE class
marker_sizes = []
for family in matched_df['TE_family']:
    if 'SINE' in family or 'LINE' in family or 'LTR' in family:
        marker_sizes.append(22)  # Marker size for Class II TEs
    elif 'DNA' in family or 'RC' in family:
        marker_sizes.append(40)  # Marker size for Class I TEs
    else:
        marker_sizes.append(0)  # Default marker size for other TEs

# Add scatter plot for TE locations with different colors and adjusted marker size based on TE class
te_scatter = ax.scatter(matched_df['chromosome_TE'], matched_df['TE_start'], color=colors, alpha=0.7, s=marker_sizes, zorder=2)

# Add scatter plot for gene locations with a valid color and adjusted transparency and marker size
ax.scatter(matched_df['Gene_chromosome'], matched_df['Gene_start'], label='Gene locations', color=diverging_colors[9], alpha=0.8, s=100)


# Add centromere markers
centromere_data = pd.read_csv('chrcen.txt', sep='\s+', header=None, names=['chromosome', 'cen_position'])
for _, row in centromere_data.iterrows():
    ax.plot(row['chromosome'], row['cen_position'], marker='D', markersize=5, color=colorblind_palette[7], zorder=2, alpha=1)

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

plt.margins(y=0.02)

ax.set_xlim(0.4, len(unique_chromosomes) + 0.6)

# Create custom legend
legend_elements = [te_scatter,
                   plt.Line2D([0], [0], marker='o', color='w', alpha=0.7, markerfacecolor=diverging_colors[0], markersize=8, zorder=0, label='DE Class II RNA TE'),
                   plt.Line2D([0], [0], marker='o', color='w', alpha=0.7, markerfacecolor=diverging_colors[6], markersize=11, label='DE Class I DNA TE'),
                   plt.Line2D([0], [0], marker='o', color='w', alpha=0.7, markerfacecolor=diverging_colors[9], markersize=15, label='DE Gene'),
                   plt.Line2D([0], [0], marker='D', color='w', alpha=1, markersize=6, markerfacecolor=colorblind_palette[7], markeredgewidth=0, label='Centromere')
                   ]

# Add legend with custom elements
ax.legend(handles=legend_elements, loc='upper right', fontsize=12)

# Show the plot
plt.tight_layout()
plt.xticks(fontsize=18, fontweight='bold')
plt.yticks(fontsize=18, fontweight='bold')
plt.savefig("te_genes_plot_phenogram_socstress.png", dpi=600)
plt.show()


# Step 4 of 5
# Plotting bar charts counting number of TE family/class counts
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def plot_te_class_counts(te_class_counts):
    """Plots the counts of each TE class."""
    sns.set_context("paper")
    fig, ax = plt.subplots(figsize=(12, 12))
    te_class_counts.plot(kind='bar', color="plum", ax=ax)
    ax.set_xlabel('TE Class', fontsize=20, fontweight='bold')
    ax.set_ylabel('Count', fontsize=20, fontweight='bold')
    ax.set_title('TE Class Counts', fontsize=20, fontweight='bold')  # Optionally adjust the title
    for tick in ax.get_xticklabels():
        tick.set_fontsize(16)
        tick.set_fontweight('bold')
        tick.set_rotation(90)  # Rotate the x-axis text by 45 degrees
    for tick in ax.get_yticklabels():
        tick.set_fontsize(16)
        tick.set_fontweight('bold')
    plt.tight_layout()  # Adjust layout to prevent overlapping
    plt.savefig('te_class_counts.png')  # Save the figure as PNG
    plt.show()

if __name__ == '__main__':
    # Read the CSV file.
    matched_df = pd.read_csv("matched_chromfilter.csv")

    # Count the occurrences of each TE class.
    te_class_counts = matched_df['TE_class'].value_counts()

    # Plot the counts of each TE class and save the figure.
    plot_te_class_counts(te_class_counts)

    print(te_class_counts)

# Step 5 of 5
# for plotting by family
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# set colorblind palette
colorblind_palette = sns.color_palette("colorblind")

def plot_te_class_counts(te_class_counts):
    """Plots the counts of each TE class."""
    sns.set_context("paper")
    fig, ax = plt.subplots(figsize=(12, 12))
    colors = []
    for te_family in te_class_counts.index:
        if ' ' in te_family in te_family:
            colors.append(colorblind_palette[2])  # Blue
        elif 'DNA' in te_family or 'Helitron' in te_family or 'hAT' in te_family\
                or 'Dada' in te_family or 'CMC-Enspm' in te_family or 'Maverick' in te_family\
                or 'PIF-Harbinger' in te_family or 'TcMar-Tc1' in te_family\
                or 'P_' in te_family: # all DNA TE P is P_ because it picks up RNA TE pao otherwise
            colors.append(colorblind_palette[3])  # Red
        else:
            colors.append(colorblind_palette[0])  # Default color
    te_class_counts.plot(kind='bar', color=colors, ax=ax)
    ax.set_xlabel('TE superfamily', fontsize=20, fontweight='bold')
    ax.set_ylabel('Count', fontsize=20, fontweight='bold')
    ax.set_title('TE superfamily Counts', fontsize=20, fontweight='bold')  # Optionally adjust the title
    for tick in ax.get_xticklabels():
        tick.set_fontsize(16)
        tick.set_fontweight('bold')
        tick.set_rotation(90)  # Rotate the x-axis text by 45 degrees
    for tick in ax.get_yticklabels():
        tick.set_fontsize(16)
        tick.set_fontweight('bold')
    plt.tight_layout()  # Adjust layout to prevent overlapping
    plt.savefig('te_class_counts.png')  # Save the figure as PNG
    plt.show()

if __name__ == '__main__':
    # Read the CSV file.
    matched_df = pd.read_csv("matched_chromfilter.csv")

    # Count the occurrences of each TE class.
    te_class_counts = matched_df['TE_class'].value_counts()

    # Plot the counts of each TE class and save the figure.
    plot_te_class_counts(te_class_counts)

    print(te_class_counts)


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

