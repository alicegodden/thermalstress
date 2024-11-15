# Title : Plotting binomial results
# Author: Dr. Alice M. Godden

import pandas as pd
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv("enrichment_test_summary_ov.csv")

# Extract relevant data for visualization
plot_data = data[data["Metric"].isin(["Observed DEGs Targeted", "Expected DEGs Targeted"])]
p_value = data[data["Metric"] == "P-Value"]["Value"].values[0]  # Extract the p-value

# Plot
plt.figure(figsize=(8, 6))
bars = plt.bar(plot_data["Metric"], plot_data["Value"], color=["#1b0c41", "#781c6d"])
# Deep Purple: #1b0c41 and # Dark Pink: #781c6d for ovaries
# for testes color=["#a52c60", "#4a0c6b"]

# Add a significance bar
x1, x2 = 0, 1  # Indices of the bars to compare
y, h, col = max(plot_data["Value"]) + 1, 0.5, "black"  # Height and offset for the significance bar
plt.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, color=col)

# Determine the significance label
if p_value < 0.001:
    significance_label = "*** (p < 0.001)"
elif p_value < 0.01:
    significance_label = "** (p < 0.01)"
elif p_value < 0.05:
    significance_label = "* (p < 0.05)"
else:
    significance_label = f"n.s. (p = {p_value:.4f})"

# Add the significance label
plt.text((x1 + x2) * 0.5, y + h, significance_label, ha="center", va="bottom", color=col, fontsize=12, fontweight="bold")

# Add title and labels
plt.title("Comparison of Observed and Expected DEGs Targeted- Ovaries", fontsize=14, fontweight="bold")
plt.ylabel("Number of DEGs", fontsize=12, fontweight="bold")
plt.xticks(fontsize=10, fontweight="bold")

# Show and save the plot
plt.tight_layout()
plt.savefig("enrichment_test_summary_with_significance_ovaries.png", dpi=300)
plt.show()
