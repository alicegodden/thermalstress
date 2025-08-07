# Title: Heatmap plotting of feature loci of TE insertion mutations from Retroseq
# Author: Dr. Alice M. Godden

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from scipy.stats import fisher_exact
import statsmodels.stats.multitest as smm
import numpy as np

# Set the path to your CSV files
csv_files_path = '*.csv'

# Custom sample name mapping
sample_name_mapping = {
    "feature_types_FC.csv": "Ovaries Control",
    "feature_types_FT.csv": "Ovaries Temperature",
    "feature_types_MC.csv": "Testes Control",
    "feature_types_MT.csv": "Testes Temperature"
}

# Get all CSV files
csv_files = glob.glob(csv_files_path)

if not csv_files:
    print("No CSV files found. Check the directory and path.")
    exit()
else:
    print(f"Found {len(csv_files)} CSV files.")

# Load and count features
feature_counts = {}
for csv_file in csv_files:
    base_filename = os.path.basename(csv_file)
    sample_name = sample_name_mapping.get(base_filename, base_filename)
    df = pd.read_csv(csv_file)
    if 'Feature' in df.columns:
        print(f"Processing {csv_file} ({len(df)} rows)")
        counts = df['Feature'].value_counts()
        feature_counts[sample_name] = counts
    else:
        print(f"Warning: 'Feature' column not found in {csv_file}")

# Create DataFrame of feature counts
feature_df = pd.DataFrame(feature_counts).fillna(0).T
feature_df = feature_df.reindex(['Ovaries Control', 'Ovaries Temperature', 'Testes Control', 'Testes Temperature'])

# --- Define the comparisons to be made ---
comparisons = {
    'Ovaries T vs. C': ('Ovaries Temperature', 'Ovaries Control'),
    'Testes T vs. C': ('Testes Temperature', 'Testes Control')
}

all_enrichment_results = []

for comparison_name, (sample1, sample2) in comparisons.items():
    print(f"\nPerforming analysis for: {comparison_name}")

    # Get data for the two samples in the comparison
    comp_df = feature_df.loc[[sample1, sample2]].copy()

    # === Enrichment Analysis for Pairwise Comparison ===
    enrichment_results = []

    sample1_total = comp_df.loc[sample1].sum()
    sample2_total = comp_df.loc[sample2].sum()

    for feature in comp_df.columns:
        a = comp_df.loc[sample1, feature]
        c = comp_df.loc[sample2, feature]
        b = sample1_total - a
        d = sample2_total - c

        contingency_table = [[a, b], [c, d]]

        try:
            oddsratio, p_value = fisher_exact(contingency_table, alternative='greater')
        except:
            oddsratio, p_value = np.nan, 1.0

        enrichment_results.append({
            'Sample': comparison_name,
            'Feature': feature,
            'OddsRatio': oddsratio,
            'PValue': p_value,
            'Count': a
        })

    enrichment_df = pd.DataFrame(enrichment_results)

    # Correct for multiple testing
    adj_pvalues = smm.multipletests(enrichment_df['PValue'], method='fdr_bh')[1]
    enrichment_df['AdjPValue'] = adj_pvalues

    all_enrichment_results.append(enrichment_df)

# Concatenate all results into a single DataFrame
final_enrichment_df = pd.concat(all_enrichment_results, ignore_index=True)

# Pivot to matrices for plotting
odds_ratio_matrix = final_enrichment_df.pivot(index='Sample', columns='Feature', values='OddsRatio')
log2_or_matrix = np.log2(odds_ratio_matrix.replace(0, np.nan))
log2_or_cleaned = log2_or_matrix.replace([np.inf, -np.inf], np.nan).fillna(0)

# Sort features by total abundance from the original dataframe
feature_totals = feature_df.sum()
sorted_features = feature_totals.sort_values().index
log2_or_cleaned = log2_or_cleaned.reindex(columns=sorted_features)

# Significance mask
sig_mask_df = final_enrichment_df.pivot(index='Sample', columns='Feature', values='AdjPValue')
sig_mask = sig_mask_df < 0.05
sig_mask = sig_mask.reindex(columns=sorted_features)

# === Plot a single Heatmap of log2(odds ratios) ===
plt.figure(figsize=(12, 4))
ax = sns.heatmap(log2_or_cleaned, cmap='vlag', annot=False, fmt=".2f",
                 cbar_kws={'label': 'log2(Odds Ratio)'}, annot_kws={"fontweight": "bold"})

cbar = ax.collections[0].colorbar
cbar.ax.yaxis.label.set_fontweight('bold')
cbar.ax.yaxis.label.set_fontsize(12)

for t in cbar.ax.yaxis.get_ticklabels():
    t.set_fontweight('bold')
    t.set_fontsize(12)

# Asterisks for significant enrichments
for y in range(log2_or_cleaned.shape[0]):
    for x in range(log2_or_cleaned.shape[1]):
        if sig_mask.iloc[y, x]:
            ax.text(x + 0.5, y + 0.5, "*", color='white', fontsize=28,
                    ha='center', va='center', fontweight='bold')

# Bold formatting
ax.set_xticklabels(ax.get_xticklabels(), fontweight='bold', fontsize=12)
ax.set_yticklabels(ax.get_yticklabels(), fontweight='bold', fontsize=12)
plt.xlabel('Feature', fontweight='bold', fontsize=12)
plt.ylabel('Comparison', fontweight='bold', fontsize=12)
plt.title("Retroseq TE Insertion Loci Features - Pairwise Enrichment Analysis", fontweight='bold', fontsize=16)

plt.tight_layout()
plt.savefig('pairwise_feature_enrichment_heatmap_log2OR.png', dpi=600)
plt.show()
