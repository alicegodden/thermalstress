# Title: TE class plots
# Author : Dr. Alice M. Godden

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Data
te_classes = ['LINE', 'LTR', 'DNA', 'SINE']
counts = [152, 4, 13, 10]

# Create plot
fig, ax = plt.subplots()

# Plot data
bars = ax.bar(te_classes, counts, color=sns.color_palette("flare", len(te_classes)))

# Set x-axis labels with bold text and rotation for readability
ax.set_xticklabels(te_classes, rotation=45, ha='right', fontweight='bold')
plt.xticks(fontweight='bold')
plt.yticks(fontweight='bold')
ax.set_xlabel("TE Class", fontweight='bold')
ax.set_ylabel("TE Class Count", fontweight='bold')
plt.title("ASM1731132v1 Sig DE TEs from Telescope by class", fontweight='bold')

# Add a legend (optional, as colors are self-explanatory with labels)
#ax.legend(bars, te_classes)

# Adjust spacing for better layout
plt.tight_layout()

# Show plot
plt.savefig('teclass_telescope_newgen.png', dpi=600)
plt.show
