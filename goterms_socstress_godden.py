import matplotlib.pyplot as plt
import numpy as np

# Data
enrichment_fdr = [0.0101080027713371,   0.00234630120792586,   0.00321341441935638,   0.0049014277017043,
                  0.00312859905463629,  0.00248485734075598,   0.00234630120792586,   0.0132748909215802,
                  0.00248485734075598,  0.00312859905463629,   0.0322838815380815,    0.0331196678119323,
                  0.000574997250239636, 0.0132748909215802]

n_genes = [3, 11, 10, 10, 14, 15, 16, 12, 16, 16, 12, 12, 24, 23]

# fold enrichment
go_terms = [26.6983240223464,   5.6477223893425,   5.33966480446927,  5.03741962685781,  3.93448985592473,
            3.88810544014753,   3.84840706628416,  3.72534753800182,  3.71454942919602,  3.53035689551688,
            3.37241987650691,   3.3372905027933,   3.28594757198109,  2.43675179569034]

# go term description
descriptions = ['Reg. of insulin-like growth factor receptor signaling pathway ',
                'Cellular component assembly involved in morphogenesis ',   'Myofibril assembly ',
                'Striated muscle cell development ',    'Muscle cell development ',
                'Striated muscle cell differentiation ',    'Striated muscle tissue development ',
                'Actomyosin structure organization ',   'Muscle tissue development ',  'Muscle cell differentiation ',
                'Muscle organ development ',    'Heart morphogenesis ',    'Muscle structure development ',
                'Heart development '

]

# Convert the descriptions list to a NumPy array
descriptions = np.array(descriptions)

# Calculate the bubble sizes
bubble_sizes = np.array(n_genes)  # set to scale to n_genes

# Create a dot plot
plt.figure(figsize=(10, 6))

# Reverse the y-axis to have each term in a row
plt.gca().invert_yaxis()

# Scatter plot with x and y data
plt.scatter(go_terms, np.arange(len(descriptions)), s=bubble_sizes, c='blue', alpha=1.0)

# Add labels and title
plt.xlabel('Fold Enrichment', fontweight='bold', fontsize=14)
plt.ylabel('GO Biological Processes', fontweight='bold', fontsize=14)
plt.title('RNA-seq: GO Biological Processes', fontweight='bold', fontsize=16, loc='center')

# Set x-ticks and labels to go_terms
plt.xticks(np.arange(0, 30, 5))
plt.yticks(np.arange(len(descriptions)), descriptions, fontweight='bold')

# Rotate the x-ticks by 90 degrees
# plt.xticks(rotation=90)

# Show the plot with extended margins
plt.margins(0.05, 0.1)

# Create a legend with bubble sizes
sizes = [5, 10, 20, 25]  # Adjust the sizes as needed
legend_labels = [plt.scatter([], [], s=size, c='blue', label=f'{size} genes') for size in sizes]

# Display the legend with the plot
legend = plt.legend(handles=legend_labels, loc='upper left', bbox_to_anchor=(1, 1), title='NGenes')
legend.get_title().set_fontweight('bold')

# Show the plot
plt.grid(True, alpha=0.25)
plt.tight_layout()
plt.savefig('dotplot_rnaseq_GO_BP_.png', dpi=600)
plt.show()

