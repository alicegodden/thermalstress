# Title: Retroseq loci in GRCz11 genome features
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
csv_files_path = '*.csv'  # Update this path if necessary

# Custom sample name mapping based on filename without extension
sample_name_mapping = {
    'OC_features.csv': 'OC',
    'OT_features.csv': 'OT',
    'TC_features.csv': 'TC',
    'TT_features.csv': 'TT'
}

# Get a list of all CSV files
csv_files = glob.glob(csv_files_path)

if not csv_files:
    print("No CSV files found. Check the directory and path.")
else:
    print(f"Found {len(csv_files)} CSV files.")

feature_counts = {}

for csv_file in csv_files:
    base_filename = os.path.basename(csv_file)
    sample_name = sample_name_mapping.get(base_filename, base_filename)
    df = pd.read_csv(csv_file)

    if 'Feature' in df.columns:
        print(f"Processing file: {csv_file} with {len(df)} rows.")
        counts = df['Feature'].value_counts()
        print(f"Feature counts for {sample_name}:\n{counts}")
        feature_counts[sample_name] = counts
    else:
        print(f"Warning: 'Feature' column not found in {csv_file}.")

if feature_counts:
    feature_df = pd.DataFrame(feature_counts).fillna(0)
    print("Feature DataFrame before transposing:\n", feature_df)
    feature_df = feature_df.T
    print("Feature DataFrame after transposing:\n", feature_df)
    feature_df = feature_df.reindex(['F0', 'F1SAT', 'F1LAT', 'F3SAT', 'F3LAT'])
    print("Feature DataFrame after reindexing:\n", feature_df)

    feature_totals = feature_df.sum(axis=0)
    sorted_features = feature_totals.sort_values().index
    filtered_feature_df = feature_df[sorted_features]
    print("Filtered Feature DataFrame for plotting:\n", filtered_feature_df)

    # === Enrichment Analysis ===
    all_feature_totals = feature_df.sum()
    enrichment_results = []

    for sample in feature_df.index:
        sample_total = feature_df.loc[sample].sum()
        for feature in feature_df.columns:
            a = feature_df.loc[sample, feature]
            b = sample_total - a
            c = all_feature_totals[feature] - a
            d = all_feature_totals.sum() - (a + b + c)
            contingency_table = [[a, b], [c, d]]

            # Fisher’s exact test (enrichment)
            try:
                oddsratio, p_value = fisher_exact(contingency_table, alternative='greater')
            except:
                oddsratio, p_value = np.nan, 1.0

            enrichment_results.append({
                'Sample': sample,
                'Feature': feature,
                'OddsRatio': oddsratio,
                'PValue': p_value,
                'Count': a
            })

    enrichment_df = pd.DataFrame(enrichment_results)
    enrichment_df['AdjPValue'] = smm.multipletests(enrichment_df['PValue'], method='fdr_bh')[1]

    # Pivot for heatmap (log2(odds ratio))
    enrichment_pivot = enrichment_df.pivot(index='Sample', columns='Feature', values='OddsRatio')
    enrichment_log2 = np.log2(enrichment_pivot.replace(0, np.nan))

    # Mask for significance
    sig_mask = enrichment_df.pivot(index='Sample', columns='Feature', values='AdjPValue') < 0.05

    # === Enrichment Heatmap ===
    plt.figure(figsize=(10, 6))
    ax = sns.heatmap(enrichment_log2, cmap='viridis', annot=True, fmt=".2f",
                     mask=enrichment_log2.isna(), cbar_kws={'label': 'log2(Odds Ratio)'})

    # Annotate significant cells
    for y in range(enrichment_log2.shape[0]):
        for x in range(enrichment_log2.shape[1]):
            if sig_mask.iloc[y, x]:
                ax.text(x + 0.5, y + 0.5, f"*", color='white', fontsize=14,
                        ha='center', va='center', weight='bold')

    plt.title('Feature Enrichment in Samples (Fisher\'s Exact Test)', fontweight='bold')
    plt.xlabel('Feature', fontweight='bold')
    plt.ylabel('Sample', fontweight='bold')
    plt.tight_layout()
    plt.savefig('feature_enrichment_heatmap_FISH_HAPLOID.png', dpi=600)
    plt.show()

    # === Stacked Bar Plot ===
    rocket_palette = sns.color_palette("rocket", n_colors=len(sorted_features))
    alpha = 0.5
    rocket_palette_with_alpha = [(r, g, b, alpha) for r, g, b in rocket_palette]

    plt.figure(figsize=(12, 6))
    filtered_feature_df.plot(kind='bar', stacked=True, figsize=(12, 6),
                             color=rocket_palette_with_alpha, edgecolor='black')

    plt.yscale('log')
    plt.title('Danio- Abundance of TE Insertion Feature Locations', fontweight='bold')
    plt.ylabel('Count (Log Scale)', fontweight='bold')
    plt.xlabel('Samples', fontweight='bold')
    plt.xticks(ticks=range(len(filtered_feature_df.index)), labels=filtered_feature_df.index, rotation=0,
               fontweight='bold')
    plt.yticks(fontweight='bold')
    plt.legend(title='Feature', labels=sorted_features, bbox_to_anchor=(1.05, 1),
               loc='upper left', title_fontsize='medium', fontsize='medium', frameon=True)
    plt.tight_layout()
    plt.savefig('retroseq_annotation_bar_chart_FISH_HAPLOID_hires.png', dpi=1200)
    plt.show()

else:
    print("No valid feature data found to plot.")
