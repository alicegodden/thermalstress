import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# Step 1: Define input files and conditions
csv_files = [
    "feature_types_FC.csv", "feature_types_FT.csv",
    "feature_types_MC.csv", "feature_types_MT.csv"
]
condition_mapping = {
    "feature_types_FC.csv": "Ovaries Control",
    "feature_types_FT.csv": "Ovaries Temperature",
    "feature_types_MC.csv": "Testes Control",
    "feature_types_MT.csv": "Testes Temperature"
}

# Combine counts from all samples
all_counts = []
for file in csv_files:
    df = pd.read_csv(file)
    condition = condition_mapping[file]
    counts = df['Feature'].value_counts().reset_index()
    counts.columns = ['Feature', 'Count']
    counts['Condition'] = condition
    all_counts.append(counts)

# Aggregate into a single DataFrame
observed_counts = pd.concat(all_counts)
observed_counts = observed_counts.pivot_table(index='Feature', columns='Condition', values='Count', aggfunc='sum').fillna(0)

# Initialize results list for Fisher's Exact Test
results = []

# Step 2: Perform Fisher's Exact Test for each feature in ovaries and testes
for feature in observed_counts.index:
    for group in ["Ovaries", "Testes"]:
        # Extract control and temperature counts
        control_col = f"{group} Control"
        temp_col = f"{group} Temperature"

        observed_control = observed_counts.loc[feature, control_col]
        observed_temperature = observed_counts.loc[feature, temp_col]

        # Total counts for the group
        total_control = observed_counts[control_col].sum()
        total_temperature = observed_counts[temp_col].sum()

        # Construct 2x2 contingency table
        contingency_table = [
            [observed_control, observed_temperature],
            [total_control - observed_control, total_temperature - observed_temperature]
        ]

        # Perform Fisher's Exact Test
        odds_ratio, p_value = fisher_exact(contingency_table)

        # Append the result
        results.append({
            "Feature": feature,
            "Group": group,
            "Observed_Control": observed_control,
            "Observed_Temperature": observed_temperature,
            "P_Value": p_value
        })

# Convert results to a DataFrame
results_df = pd.DataFrame(results)

# Step 3: Correct p-values for multiple testing (FDR)
results_df['P_Value_Corrected'] = multipletests(results_df['P_Value'], method='fdr_bh')[1]

# Save the results
results_df.to_csv("fisher_results_full.csv", index=False)

# Filter significant features
significant_features = results_df[results_df['P_Value_Corrected'] < 0.05]
significant_features.to_csv("fisher_results_significant.csv", index=False)

# Step 4: Plot the results
plt.figure(figsize=(15, 8))

# Extract conditions and features
conditions = ["Ovaries Control", "Ovaries Temperature", "Testes Control", "Testes Temperature"]
features = observed_counts.index

# Bar width and spacing
bar_width = 0.9
bar_spacing = 0.3
x_positions = np.arange(len(features)) * len(conditions) * bar_spacing

# Apply a color palette
palette = sns.color_palette("rocket", n_colors=len(conditions))

# Plot each condition
for j, condition in enumerate(conditions):
    bar_positions = x_positions + j * bar_spacing
    plt.bar(
        bar_positions, observed_counts[condition], width=bar_width / len(conditions),
        label=condition, color=palette[j], edgecolor='black'
    )

# Add significance stars
for _, row in significant_features.iterrows():
    feature = row["Feature"]
    group = row["Group"]
    condition = f"{group} Temperature"  # Always for temperature condition

    # Get bar position and height
    feature_index = list(features).index(feature)
    bar_index = conditions.index(condition)
    bar_pos = x_positions[feature_index] + bar_index * bar_spacing
    bar_height = observed_counts.loc[feature, condition]

    # Add a star above the bar
    plt.text(
        bar_pos, bar_height * 1.05, '*',
        ha='center', fontsize=14, fontweight='bold', color='black'
    )

# Set y-axis to log scale
plt.yscale('log')

# Add labels and title
plt.xlabel('Feature', fontweight='bold', fontsize=18)
plt.ylabel('Count (Log Scale)', fontweight='bold', fontsize=18)
plt.yticks(fontweight="bold", fontsize=14)
plt.title('TE Insertion Enrichment by Feature Type', fontweight='bold', fontsize=20)

# Add x-ticks
xtick_positions = x_positions + ((len(conditions) - 1) * bar_spacing) / 2
plt.xticks(ticks=xtick_positions, labels=features, rotation=45, ha='right', fontsize=14, fontweight='bold')

# Add vertical boundary lines
for i in range(len(features) - 1):
    plt.axvline(
        (xtick_positions[i] + xtick_positions[i + 1]) / 2, color='grey', linestyle='--', linewidth=0.8, alpha=0.7
    )

# Add legend outside the plot
plt.legend(
    title='Treatment', fontsize=14, title_fontsize=16, bbox_to_anchor=(1.05, 1), loc='upper left',
    prop={'weight': 'bold'}  # Bold legend text
)

# Adjust layout and save
plt.tight_layout()
plt.savefig("fisher_bar_chart_with_significance.png", dpi=300)
plt.show()
