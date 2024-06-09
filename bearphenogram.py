# Title: Bear phenogram for plotting sig DE TEs from telescope outputs
# Author : Dr. Alice M. Godden

# Notes
# Used grep on command line along with the TE.rmsk.gtf annotation file to find genomic co-ordinates of sig DE TEs, saved this as a list .csv file.

#Import libraries 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read chromosome end data
chrom_end_df = pd.read_csv("chrom_end_bear.txt", header=None, names=["Chromosome", "Length"], delim_whitespace=True)

# Read points data
points_df = pd.read_csv("output_matches_TEs_trim.txt", names=["Chromosome", "Start", "End", "Name"],
                        delim_whitespace=True)

# Get unique chromosomes with hits
chromosomes_with_hits = points_df["Chromosome"].unique()

# Filter chromosome end data for chromosomes with hits and sort by chromosome number
chrom_end_df_filtered = chrom_end_df[chrom_end_df["Chromosome"].isin(chromosomes_with_hits)]
chrom_end_df_filtered = chrom_end_df_filtered.sort_values(by="Chromosome")

# Define colors for different element types
try:
    flare_palette = sns.color_palette("rocket")
    element_colors = {"DNA": flare_palette[0],
                      "LINE": flare_palette[2],
                      "SINE": flare_palette[4],
                      "LTR": flare_palette[5]}
except IndexError:
    print("Error: Flare palette does not contain enough colors.")
    element_colors = {"DNA": "blue",
                      "LINE": "green",
                      "SINE": "red",
                      "LTR": "purple"}

# Plot chromosomes with hits and corresponding points
plt.figure(figsize=(10, 6))
legend_handles = []
legend_labels = set()  # Set to store legend labels
for index, row in chrom_end_df_filtered.iterrows():
    chrom = row["Chromosome"]
    length = row["Length"]
    plt.plot([0, length], [chrom, chrom], color="black")
    points_on_chrom = points_df[points_df["Chromosome"] == chrom]
    for _, point_row in points_on_chrom.iterrows():
        start = point_row["Start"]
        end = point_row["End"]
        element_type = point_row["Name"].split("#")[1].split("/")[0]  # Extract element type from the Name column
        color = element_colors.get(element_type, "gray")  # Use gray color if element type is not found in the dictionary
        plt.scatter([start, end], [chrom, chrom], color=color, alpha=0.6)
        # Add to legend only if it hasn't been added before
        if element_type not in legend_labels:
            legend_handles.append(plt.scatter([], [], color=color, label=element_type, alpha = 0.6))
            legend_labels.add(element_type)

plt.xlabel("Genomic Position", fontweight='bold')
plt.ylabel("Chromosome/Scaffold", fontweight='bold')
plt.title("Genomic Positions of significantly differentially expressed TEs", fontweight='bold')
plt.yticks(chrom_end_df_filtered["Chromosome"], chrom_end_df_filtered["Chromosome"])
plt.yticks(fontweight='bold')
plt.xticks(fontweight='bold')
plt.legend(handles=legend_handles, title="TE Class")
plt.gca().invert_yaxis()  # Invert y-axis to plot chromosomes from top to bottom
plt.grid(False)  # Remove gridlines
plt.show()
