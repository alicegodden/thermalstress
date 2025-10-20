# Title: Plotting Retroseq results processed with Ensembl VEP
# Author: Dr. Alice M. Godden

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import fisher_exact  # Added for statistical testing

# --- Configuration ---
# NOTE: Replace 'vep_annotated_MEIs_fixed.tsv' with the actual path to your VEP output file
VEP_INPUT_FILE = 'TT.vep.ens.txt'
PLOT_SUMMARY_OUTPUT = 'vep_consequence_summary_TT_gq100_fl8.png'  # Filename retained from user request
PLOT_HEATMAP_OUTPUT = 'vep_gene_consequence_heatmap_TT_gq100_fl8.png'

# Scientific order for IMPACT for proper plotting (descending severity)
IMPACT_ORDER = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']


# --- 1. Data Loading and Initial Cleaning ---

def load_and_clean_data(file_path):
    """
    Loads VEP output, performs essential column mapping, and extracts the most severe consequence.
    Returns the full processed DataFrame, which is then filtered downstream for specific plots.
    """
    print(f"Loading VEP output from: {file_path}")

    try:
        # Read the tab-separated file
        df = pd.read_csv(file_path, sep='\t', low_memory=False)
    except FileNotFoundError:
        print(f"ERROR: Input file not found at '{file_path}'. Please check the path.")
        return None

    # --- CRITICAL: Map VEP's column names (SYMBOL, IMPACT, Consequence) ---
    column_map = {}
    if 'Consequence' in df.columns:
        column_map['Consequence'] = 'Consequence_Raw'
    if 'SYMBOL' in df.columns:
        column_map['SYMBOL'] = 'Gene_Symbol'
    if 'IMPACT' in df.columns:
        column_map['IMPACT'] = 'Impact_Level'

    if not all(col in df.columns for col in ['Consequence', 'SYMBOL', 'IMPACT']):
        missing_cols = [col for col in ['Consequence', 'SYMBOL', 'IMPACT'] if col not in df.columns]
        print(f"ERROR: Required columns are missing: {', '.join(missing_cols)}. Check your header row.")
        return None

    df = df.rename(columns=column_map)

    # --- VEP lists all consequences. We take the first one (MOST SEVERE).
    df['Most_Severe_Consequence'] = df['Consequence_Raw'].apply(
        lambda x: str(x).split('&')[0].split(',')[0] if pd.notna(x) else 'Unknown'
    )

    # Only filter by valid IMPACT levels here. We keep all rows, even those without a gene.
    df_processed = df[
        df['Impact_Level'].isin(IMPACT_ORDER)  # Ensure Impact is one of the four valid levels
    ].copy()  # Use .copy() to avoid SettingWithCopyWarning

    print(f"Loaded and processed {len(df_processed)} variants with valid impact levels.")
    return df_processed


# --- 2. Plotting Functions ---

def plot_consequence_summary(df, output_file):
    """
    Generates and saves a bar plot of ALL consequences, colored by IMPACT,
    including statistical annotation for HIGH impact enrichment.
    """

    sns.set_style("whitegrid")

    # Group by Consequence and Impact to get counts
    plot_df = df.groupby(['Most_Severe_Consequence', 'Impact_Level']).size().reset_index(name='Count')
    total_variants = plot_df['Count'].sum()

    # --- STATISTICAL ANALYSIS: Enrichment of HIGH Impact ---

    P_THRESHOLD = 0.05
    OR_THRESHOLD = 1.0
    enrichment_results = {}

    # 1. Calculate overall counts for the reference background
    all_high_counts = plot_df[plot_df['Impact_Level'] == 'HIGH']['Count'].sum()
    all_non_high_counts = total_variants - all_high_counts

    for consequence in plot_df['Most_Severe_Consequence'].unique():
        con_df = plot_df[plot_df['Most_Severe_Consequence'] == consequence]

        # a: HIGH impact count in current consequence (Test Group HIGH)
        a = con_df[con_df['Impact_Level'] == 'HIGH']['Count'].sum()
        # b: Non-HIGH impact count in current consequence (Test Group Non-HIGH)
        b = con_df[con_df['Impact_Level'] != 'HIGH']['Count'].sum()

        # c: HIGH impact count in all OTHER consequences (Background HIGH)
        c = all_high_counts - a
        # d: Non-HIGH impact count in all OTHER consequences (Background Non-HIGH)
        d = all_non_high_counts - b

        table = [[a, b], [c, d]]

        if (a + b) == 0 or (c + d) == 0:
            enrichment_results[consequence] = (1.0, 1.0, 0)
            continue

        try:
            # alternative='greater' tests for enrichment (Odds Ratio > 1)
            oddsratio, p_value = fisher_exact(table, alternative='greater')
            enrichment_results[consequence] = (oddsratio, p_value, a + b)
        except ValueError:
            enrichment_results[consequence] = (np.nan, 1.0, a + b)

    # --- END STATISTICAL ANALYSIS ---

    # Calculate Total Count per Consequence (for sorting)
    consequence_totals = plot_df.groupby('Most_Severe_Consequence')['Count'].sum().sort_values(
        ascending=False).index.tolist()

    # 1. Sort the bars on the X-axis by IMPACT, then by total count
    plot_df['Impact_Level'] = pd.Categorical(plot_df['Impact_Level'], categories=IMPACT_ORDER, ordered=True)

    # Create a sort index based on total counts for visual hierarchy
    consequence_sorter = {con: i for i, con in enumerate(consequence_totals)}
    plot_df['Consequence_Sort_Index'] = plot_df['Most_Severe_Consequence'].apply(
        lambda x: consequence_sorter.get(x, len(consequence_totals)))

    # Final sort for plotting: severe impact consequences first, then by total count
    plot_df.sort_values(by=['Impact_Level', 'Consequence_Sort_Index'], ascending=[True, True], inplace=True)

    # Define the final X order for the plot
    X_ORDER = plot_df['Most_Severe_Consequence'].unique()

    plt.figure(figsize=(10, 8))
    ax = sns.barplot(
        x='Most_Severe_Consequence',
        y='Count',
        data=plot_df,
        hue='Impact_Level',
        hue_order=IMPACT_ORDER,
        order=X_ORDER,
        palette='rocket',
        dodge=False,
        width=0.7,
        ax=plt.gca()
    )

    # --- STATISTICAL ANNOTATION: Add asterisks for significant enrichment ---
    max_heights = plot_df.groupby('Most_Severe_Consequence')['Count'].sum().to_dict()
    x_pos = dict(zip(X_ORDER, ax.get_xticks()))

    for consequence in X_ORDER:
        odds_ratio, p_value, total_count = enrichment_results.get(consequence, (1.0, 1.0, 0))

        if p_value < P_THRESHOLD and odds_ratio > OR_THRESHOLD:
            if p_value < 0.001:
                sig_text = "***"
            elif p_value < 0.01:
                sig_text = "**"
            else:
                sig_text = "*"

            x = x_pos[consequence]
            y = max_heights[consequence]
            buffer = max(5, y * 0.02)

            ax.text(x, y + buffer, sig_text,
                    ha='center', va='bottom',
                    fontsize=16, fontweight='bold', color='black')

    # --- Customization for Publication Quality ---
    plt.title(
        'Testes Temp: Distribution of All Variant Consequences, Colored by IMPACT (N={})\n'
        'Enrichment of HIGH impact variants (Fisher\'s Exact Test, one-sided)'.format(
            total_variants),
        fontsize=14, fontweight='bold', pad=20
    )
    plt.xlabel('Variant Consequence', fontsize=14, fontweight='bold')
    plt.ylabel('Number of Variants', fontsize=14, fontweight='bold')
    plt.yscale('log')

    plt.xticks(rotation=45, ha='right', fontsize=14, fontweight='bold')
    plt.yticks(rotation=0, ha='right', fontsize=14, fontweight='bold')

    legend = ax.legend(title='Impact Level', bbox_to_anchor=(1.01, 1), loc='upper left')
    legend.get_title().set_fontweight('bold')

    for text in legend.get_texts():
        text.set_fontweight('bold')
        text.set_fontsize(12)

    sns.despine(left=False, bottom=False)
    plt.tight_layout(rect=[0, 0.03, 0.95, 0.9])
    plt.savefig(output_file, dpi=300)
    print(f"\nSummary Bar Plot saved to: '{output_file}'")
    plt.close()


def plot_consequence_heatmap(df, output_file, top_n=20):
    """
    Generates and saves a heatmap showing consequence distribution across top genes,
    with statistical enrichment annotations for each cell.
    """

    # Filter the DataFrame to include ONLY variants with valid gene symbols
    df_heatmap_full = df[
        df['Gene_Symbol'].notna() &
        (df['Gene_Symbol'] != '') &
        (df['Gene_Symbol'] != '-')
        ].copy()

    print(
        f"Generating heatmap for the Top {top_n} most mutated genes (from {len(df_heatmap_full)} gene-associated variants)...")

    if df_heatmap_full.empty:
        print("Not enough gene-associated data to create a heatmap. Skipping.")
        return

    # 1. Identify the top N most frequently mutated genes
    top_genes = df_heatmap_full['Gene_Symbol'].value_counts().nlargest(top_n).index.tolist()
    df_top_genes = df_heatmap_full[df_heatmap_full['Gene_Symbol'].isin(top_genes)]

    if df_top_genes.empty:
        print("Not enough data to create a heatmap with valid genes. Skipping.")
        return

    # 2. Create the contingency table (pivot table of counts) for plotting
    heatmap_data = pd.crosstab(
        df_top_genes['Gene_Symbol'],
        df_top_genes['Most_Severe_Consequence']
    )

    # 3. Determine the desired column order for the heatmap (sort by total counts for clarity)
    consequence_totals = heatmap_data.sum(axis=0).sort_values(ascending=False).index.tolist()
    heatmap_data = heatmap_data.reindex(columns=consequence_totals, fill_value=0)  # Reindex and fill missing with 0

    # --- STATISTICAL ANALYSIS: Enrichment of Consequence Y in Gene X ---

    grand_total = df_heatmap_full.shape[0]  # Total number of gene-associated variants
    p_value_matrix = pd.DataFrame(1.0, index=heatmap_data.index, columns=heatmap_data.columns)

    for gene in heatmap_data.index:
        for consequence in heatmap_data.columns:
            # a: Count of Consequence Y in Gene X (cell value)
            a = heatmap_data.loc[gene, consequence]

            # b: Count of All Other Consequences in Gene X (Row Total - a)
            b = df_top_genes[df_top_genes['Gene_Symbol'] == gene].shape[0] - a

            # c: Count of Consequence Y in All Other Genes (Column Total - a)
            c = df_heatmap_full[df_heatmap_full['Most_Severe_Consequence'] == consequence].shape[0] - a

            # d: Count of All Other Consequences in All Other Genes (Grand Total - a - b - c)
            d = grand_total - a - b - c

            table = [[a, b], [c, d]]

            if (a + b) > 0 and (c + d) > 0 and (a + c) > 0 and (b + d) > 0:
                try:
                    # alternative='greater' tests for enrichment (Odds Ratio > 1)
                    oddsratio, p_value = fisher_exact(table, alternative='greater')
                    # Store p-value only if enrichment is detected (OR > 1)
                    if oddsratio > 1.0:
                        p_value_matrix.loc[gene, consequence] = p_value
                except ValueError:
                    pass

    # --- Create Plot ---
    P_THRESHOLD = 0.05
    plt.figure(figsize=(10, max(8, len(top_genes) * 0.5)))  # Dynamic height

    ax = sns.heatmap(
        heatmap_data,
        annot=False,  # Show the count numbers (Set to True for counts)
        fmt='d',  # Format as integers
        cmap='rocket_r',  # Color scheme (reverse rocket)
        linewidths=0.5,  # Lines between cells
        linecolor='black',
    )

    # --- STATISTICAL ANNOTATION on Heatmap ---
    for i, gene in enumerate(heatmap_data.index):
        for j, consequence in enumerate(heatmap_data.columns):
            p_value = p_value_matrix.loc[gene, consequence]

            if p_value < P_THRESHOLD:
                # Determine annotation text based on p-value
                if p_value < 0.001:
                    sig_text = "***"
                elif p_value < 0.01:
                    sig_text = "**"
                else:
                    sig_text = "*"

                # Place the annotation in the center of the cell
                ax.text(j + 0.5, i + 0.5, sig_text,
                        ha='center', va='center',
                        fontsize=14, fontweight='bold', color='white')

    # --- Customization for Publication Quality ---

    # FIX: Set Color Bar Label
    cbar = ax.collections[0].colorbar
    cbar.set_label('Number of Variants', weight='bold', fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    # The previous font change for yticks on the cbar often fails, simplifying here:
    cbar.ax.set_yticklabels([label.get_text() for label in cbar.ax.get_yticklabels()], weight='bold')
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=14, fontweight='bold')
    ax.set_yticklabels(ax.get_yticklabels(), ha='right', fontsize=14, fontweight='bold')

    plt.title(
        f'Testes Temp: Consequence Distribution Across Top {len(top_genes)} Mutated Genes\n'
        '(Fisher\'s Exact Test, one-sided) for the specific consequence in the gene',
        fontsize=16,
        fontweight='bold',
        pad=20
    )
    plt.ylabel('Gene Symbol', fontsize=16, fontweight='bold')
    plt.xlabel('Most Severe Consequence', fontsize=16, fontweight='bold')

    # Rotate x-axis labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=14, fontweight='bold')
    ax.set_yticklabels(ax.get_yticklabels(), ha='right', fontsize=14, fontweight='bold')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    # Save the figure
    plt.savefig(output_file, dpi=300)
    print(f"\nHeatmap saved successfully to: '{output_file}'")
    plt.close()


# --- 3. Main Execution ---
if __name__ == "__main__":

    # Load and clean returns the full DataFrame (including non-gene variants)
    processed_df = load_and_clean_data(VEP_INPUT_FILE)

    if processed_df is not None and not processed_df.empty:
        # Summary Plot: Use the full processed DataFrame
        plot_consequence_summary(processed_df, PLOT_SUMMARY_OUTPUT)

        # Heatmap Plot: The gene-filtering is now handled INSIDE this function
        # Using a default of top_n=20 as specified in your original script
        plot_consequence_heatmap(processed_df, PLOT_HEATMAP_OUTPUT, top_n=20)
    else:
        print("Analysis stopped. No valid data to plot.")
