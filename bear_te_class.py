# Title: TE class plots
# Author: Dr. Alice M. Godden

import matplotlib.pyplot as plt
import seaborn as sns

# Data
te_classes = ['LINE', 'LTR', 'DNA', 'SINE']
counts = [38, 7, 10, 8]
#countsursmar = 38, 7, 10, 8
# countsnewgen 152, 4, 13, 10

# Define color palette
flare_palette = sns.color_palette("rocket")
element_colors = {"DNA": flare_palette[0],
                  "LINE": flare_palette[2],
                  "SINE": flare_palette[4],
                  "LTR": flare_palette[5]}

# Create plot
fig, ax = plt.subplots()

# Plot data with specified colors
bars = ax.bar(te_classes, counts, color=[element_colors[te] for te in te_classes])

# Set x-axis labels with bold text and rotation for readability
ax.set_xticklabels(te_classes, rotation=45, ha='right', fontweight='bold')
ax.set_xlabel("TE Class", fontweight='bold')
ax.set_ylabel("TE Class Count", fontweight='bold')
plt.title("ASM1731132v1 Sig DE TEs from Telescope by class", fontweight='bold')

# Set the ticks fontweight
plt.xticks(fontweight='bold')
plt.yticks(fontweight='bold')

# Adjust spacing for better layout
plt.tight_layout()

# Save and show plot
plt.savefig('teclass_telescope_newgen2.png', dpi=600)
plt.show()
