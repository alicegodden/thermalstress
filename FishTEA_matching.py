# Title: FishTEA matching overlapping TEs and genes
# Subtitle: Step 2 of 2 in python FishTEA
# Author: Dr. Alice M. Godden

import pandas as pd

# Define file paths
sigTEs_file = "sigTEs.csv"
sigGenes_file = "sigGenes.csv"
log_file = "matching_log.txt"

# Read the CSV files
print("Reading CSV files...")
sigTEs_df = pd.read_csv(sigTEs_file)
sigGenes_df = pd.read_csv(sigGenes_file)

# Create an empty list to store the matched TEs and genes
matched_records = []

# Iterate over each row in sigTEs_df
print("Matching TEs with genes...")
for index_TE, TE_row in sigTEs_df.iterrows():
    TE_start = TE_row['TE_start']
    TE_end = TE_row['TE_end']
    TE_chromosome = TE_row['chromosome_TE']
    
    # Iterate over each row in sigGenes_df
    for index_gene, gene_row in sigGenes_df.iterrows():
        gene_start = gene_row['Gene_start']
        gene_end = gene_row['Gene_end']
        gene_chromosome = gene_row['Gene_chromosome']
        
        # Check if the chromosome information matches
        if TE_chromosome == gene_chromosome:
            # Check for overlap between TE and gene
            if (TE_start <= gene_end) and (gene_start <= TE_end):
                # Store the matched TE and gene information
                matched_records.append({
                    'chromosome_TE': TE_chromosome, # Include chromosome_TE in the matched records
                    'TE_name': TE_row['TE_name'],
                    'TE_class': TE_row['TE_class'],
                    'TE_family': TE_row['TE_family'],
                    'TE_start': TE_start,
                    'TE_end': TE_end,
                    'Gene_ID': gene_row['Gene_ID'],
                    'Gene_name': gene_row['Gene_name'],
                    'Gene_start': gene_start,
                    'Gene_end': gene_end,
                    'Gene_chromosome': gene_chromosome
                })
    
    # Print progress update
    print(f"Processed {index_TE + 1} out of {len(sigTEs_df)} TEs.", end="\r")

# Create a DataFrame from the matched_records list
matched_df = pd.DataFrame(matched_records)

# Write the matched dataframe to a new CSV file
matched_df.to_csv("matched.csv", index=False)

# Write log information to a log file
with open(log_file, "w") as f:
    f.write("Matching TEs with genes completed.\n")
    f.write(f"Number of TEs processed: {len(sigTEs_df)}\n")
    f.write(f"Number of matched TEs and genes: {len(matched_df)}\n")

print("\nMatching completed. Results saved to 'matched.csv'.")
