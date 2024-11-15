# Title : Binomial testing enrichment of genes targeted by miRNAs
# Author: Dr. Alice M. Godden

from scipy.stats import binomtest
import pandas as pd

# Define parameters
k = 6  # Number of DEGs targeted by DE miRNAs
n = 340  # Total number of DEGs in ovaries
p = 3 / 374  # Probability of random targeting by DE miRNAs (41 DE miRNAs out of 374 total)

# Calculate expected value
expected = n * p

# Perform the binomial test
result = binomtest(k, n, p, alternative='greater')

# Print result
print(f"P-value: {result.pvalue}")
if result.pvalue < 0.05:
    print("The result is statistically significant: DE miRNAs target DEGs more than by chance.")
else:
    print("The result is not statistically significant: No evidence for enrichment.")

# Create a summary table
summary_table = pd.DataFrame({
    "Metric": ["Observed DEGs Targeted", "Expected DEGs Targeted", "Total DEGs", "Probability of Targeting", "P-Value"],
    "Value": [k, expected, n, p, result.pvalue]
})

# Save table to a CSV and display
summary_table.to_csv("enrichment_test_summary_te.csv", index=False)
print("\nSummary Table:")
print(summary_table)

# num genes in ov background 20207
# num genes in te background 19619
# num mirnas in total is 374
# num mirnas in ov 41
# num mirnas in te 3
# num genes in list targeted by mirnas ovaries is 304
# num genes in list targeted by mirnas testes is 6
# num sig de genes in ov is 378
# num sig de genes in te is 340
