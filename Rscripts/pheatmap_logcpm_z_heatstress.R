# Title: Heatmap of normalised counts logCPM 
# Author: Dr. Alice M. Godden

# Load necessary libraries
library(DESeq2)
library(pheatmap)
library(viridis) # For the inferno color palette
library(stringr)
library(readr)

# --- Step 1: Load and process the normalized counts data ---
# Load normalized counts from a CSV file.
# The `row.names = 1` argument is crucial to use the gene IDs as row names.
counts_file <- "HS_allsamples_normalised_counts_piRNAs.csv"
normalized_counts <- read.csv(counts_file, row.names = 1)

# Calculate CPM (Counts Per Million) to normalize for library size.
# The `sweep` function divides each column by its total sum and then scales by 1e6.
cpm <- sweep(normalized_counts, 2, colSums(normalized_counts), FUN = "/") * 1e6

# Log-transform the CPM values with a pseudo-count of 1 to handle zeros.
logCPM <- log2(cpm + 1)

# Z-score normalize the data across samples for each gene.
# This step is essential for creating a heatmap where rows (genes) have comparable scales.
logCPM_z <- t(scale(t(logCPM)))

# --- Step 2: Prepare the data for the heatmap ---
# Assign the full Z-score normalized data to the heatmap_data object.
# The script no longer filters for a specific list of genes.
heatmap_data <- logCPM_z

# --- Step 3: Define annotations and labels for the heatmap ---
# These variables were previously undefined and caused the script to fail.
# I've defined them here based on the sample names from your `head(countData)` output.

# Define custom column labels for the heatmap, matching the column order of your data.
# This assumes your data has 24 columns, corresponding to 4 groups of 6 samples each.
custom_column_labels <- c(
  "OC1", "OC2", "OC3", "OC4", "OC5", "OC6",
  "OT1", "OT2", "OT3", "OT4", "OT5", "OT6",
  "TC1", "TC2", "TC3", "TC4", "TC5", "TC6",
  "TT1", "TT2", "TT3", "TT4", "TT5", "TT6"
)

# Define treatment labels for annotation, matching the order of the custom labels.
# Based on the sample names: OC (Ovaries Control), OT (Ovaries Treated),
# TC (Testes Control), TT (Testes Treated).
treatment_labels <- c(
  rep("Ovaries_Control", 6), rep("Ovaries_Temp", 6),
  rep("Testes_Control", 6), rep("Testes_Temp", 6)
)

# Create an annotation data frame. The `rownames` must match the `custom_column_labels`.
annotation_col <- data.frame(
  Treatment_Organ = factor(treatment_labels, levels = unique(treatment_labels))
)
rownames(annotation_col) <- custom_column_labels

# Define the colors for each treatment group in the annotation bar.
annotation_colors <- list(
  Treatment_Organ = c(
    "Ovaries_Control" = "#009E73",
    "Ovaries_Temp" = "#E69F00",
    "Testes_Control" = "#56B4E9",
    "Testes_Temp" = "#CC79A7"
  )
)

# --- Step 4: Plot the heatmap with the prepared data and annotations ---
# Ensure the columns in your heatmap data match the labels and annotations.
if (ncol(heatmap_data) != length(custom_column_labels)) {
  stop("The number of columns in heatmap_data does not match the number of labels.")
}

# --- Step 5: Save the heatmap plot with specific dimensions and resolution ---
# To save the plot with specific dimensions and a resolution (ppi), we use the png() function.
# Here, we set the dimensions to 10 inches by 12 inches at 300 ppi.
png(filename = "piRNAs_HS_normalised_counts.png", 
    width = 5, 
    height = 6, 
    units = "in", 
    res = 300)

pheatmap(heatmap_data, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, # Do not cluster columns; keep the original order
         show_rownames = FALSE, # Set to FALSE to avoid clutter with many genes
         show_colnames = TRUE, 
         fontface = "bold",
         scale = "none", # Data is already Z-score normalized
         color = inferno(256),
         labels_col = custom_column_labels,
         annotation_col = annotation_col, # Add the column annotations
         annotation_colors = annotation_colors, # Apply colors for each treatment group
         annotation_legend = TRUE,
         main = "Z-Score Normalized logCPM for all piRNAs")


dev.off() # Close the plotting device to finalize the file.
