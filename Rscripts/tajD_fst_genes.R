# Title : FST and Tajima's D Chromosome level analysis
# Author : Dr. Alice M. Godden

# --- IMPORTANT: Run all library() calls first! ---
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel) # For non-overlapping text labels
library(viridis) # For color palettes, though not directly used for genes in this version
# ~/Desktop/Projects/Haploid selection project/2025/fst_tajd
# --- Load data ---
# Read the Tajima's D data
file_path_tajima <- "Taj_MT_MC.Tajima.D"
tajima_data <- read.table(file_path_tajima, header = TRUE, sep = "\t")

# Read the FST data
file_path_fst <- "fst_outputMTvMC.windowed.weir.fst"
fst_data <- read.table(file_path_fst, header = TRUE, sep = "\t")

# Load gene data
# IMPORTANT: The provided gene data file 'GRCz11_genes_name_chr1_tabdelim.csv'
# is specific to chromosome 1. If your 'target_chromosome' below is not '1',
# then very few or no genes will be found and plotted from this file.
# Please ensure your gene data matches the chromosome you are plotting.
# MODIFICATION: Read in 'GRCz11_ENSEMBL_BIOMART_Chr4_genes.txt'
# MODIFICATION: Replace empty cells (read as NA) with 0 in numeric columns
gene_data_red <- read.table("chr18_genes.txt", header = TRUE,
                            sep = "\t", # Added: Specify tab as the separator
                            na.strings = c("", "NA")) %>% # Treat empty strings and "NA" as NA
  mutate(across(where(is.numeric), ~replace_na(., 0))) # Replace NAs with 0 in all numeric columns


# Filter to only chromosomes 1-25 and convert CHROM to factor for correct ordering
chromosomes_to_plot <- as.character(1:25)
tajima_data <- tajima_data %>% filter(CHROM %in% chromosomes_to_plot)
fst_data <- fst_data %>% filter(CHROM %in% chromosomes_to_plot)

tajima_data$CHROM <- factor(tajima_data$CHROM, levels = chromosomes_to_plot)
fst_data$CHROM <- factor(fst_data$CHROM, levels = chromosomes_to_plot)

# Merge Tajima's D and FST data by CHROM and BIN_START
merged_data <- full_join(
  tajima_data %>% select(CHROM, BIN_START, TajimaD),
  fst_data %>% select(CHROM, BIN_START, WEIGHTED_FST),
  by = c("CHROM", "BIN_START")
)

# --- CORRECTED SCALING FOR FST_SCALED FOR 0-POINT ALIGNMENT ---
merged_data <- merged_data %>%
  mutate(FST_SCALED = WEIGHTED_FST / 0.05) # This scales FST to the TajimaD axis with 0-alignment

# --- FST_GROUP now uses actual WEIGHTED_FST values for clear thresholding ---
merged_data <- merged_data %>%
  mutate(
    FST_GROUP = case_when(
      WEIGHTED_FST > 0.1 ~ "FST > 0.1", # Brighter orange, based on actual FST value
      TRUE ~ "FST" # Yellowish orange (FST <= 0.1)
    )
  )

# Define centromere positions based on provided data
centromere_data <- data.frame(
  CHROM = factor(1:25, levels = 1:25), # Ensure factor order matches
  CENTROMERE_POS = c(33133433, 19743539, 49120158, 25975534, 50391699,
                     34866988, 60898697, 30099230, 13283057, 9894335,
                     9507769, 27779269, 19555524, 16933009, 13943685,
                     19463555, 48225166, 24528805, 19777257, 11319545,
                     28908825, 21520944, 13631408, 15196971, 20793603)
)

# --- FILTER FOR TARGET CHROMOSOME ---
target_chromosome <- "18" # Set your target chromosome here
merged_data_filtered <- merged_data %>% filter(CHROM == target_chromosome)
centromere_data_filtered <- centromere_data %>% filter(CHROM == target_chromosome)

# --- Prepare Gene Data for Plotting ---
# Filter gene data to the target chromosome
genes_on_target_chr <- gene_data_red %>%
  filter(chrom == target_chromosome) %>% # FIX: Changed 'chromosome' to 'chrom'
  filter(!attribute %in% c("NA", "havana", "ensembl_havana", "ensembl")) %>%
  distinct(attribute, .keep_all = TRUE) # Remove duplicate gene names

# DIAGNOSTIC: Check how many genes are found for the target chromosome
message(paste("Found", nrow(genes_on_target_chr), "genes on chromosome", target_chromosome, "after initial filtering."))


# Identify the 100kb windows that contain "highlighted SNPs" (FST > 0.1)
highlighted_windows_ranges <- merged_data_filtered %>%
  filter(WEIGHTED_FST > 0.1) %>%
  mutate(window_start = BIN_START,
         window_end = BIN_START + 99999) # Assuming 100kb windows

# DIAGNOSTIC: Check how many highlighted FST windows are found
message(paste("Found", nrow(highlighted_windows_ranges), "highlighted FST windows (FST > 0.1) on chromosome", target_chromosome, "."))


# Find genes that overlap any of these highlighted FST windows
genes_to_plot_final <- data.frame()
if (nrow(highlighted_windows_ranges) > 0) {
  for (i in 1:nrow(genes_on_target_chr)) {
    gene_start <- genes_on_target_chr$start[i]
    gene_end <- genes_on_target_chr$end[i]
    
    # Check for overlap: gene starts before window_end AND gene ends after window_start
    overlaps <- highlighted_windows_ranges %>%
      filter(
        (gene_start <= window_end & gene_end >= window_start)
      )
    
    if (nrow(overlaps) > 0) {
      genes_to_plot_final <- bind_rows(genes_to_plot_final, genes_on_target_chr[i, ])
    }
  }
}
# Ensure unique genes are kept if filtering results in duplicates from multiple overlaps
genes_to_plot_final <- distinct(genes_to_plot_final, attribute, .keep_all = TRUE)

# DIAGNOSTIC: Check how many genes are finally selected for plotting
message(paste("Selected", nrow(genes_to_plot_final), "genes that overlap highlighted FST windows for plotting."))


# Define the Y-position for the gene tracks and labels
gene_track_y_position_all_genes <- -5.3 # For all genes track (lower than highlighted)
gene_track_y_position_all_genes_end <- -5.1 # To give segments a slight vertical extent

# Calculate staggered y-positions for highlighted gene labels
if (nrow(genes_to_plot_final) > 0) {
  # Base Y position for labels, ensuring they are between -5 and -2
  base_label_y_upper_limit <- -2.0 # Highest point labels should go
  min_vertical_separation <- 0.3 # Minimum vertical space between labels (e.g., for each unique Y-level)
  
  genes_to_plot_final <- genes_to_plot_final %>%
    arrange(start) %>% # Sort by start position for consistent staggering
    mutate(
      # Assign a unique "row" index for vertical staggering, forcing each label to a different y-level
      stagger_row_index = row_number(),
      # Calculate the base y for this row
      # Labels will descend from base_label_y_upper_limit, with min_vertical_separation between each
      label_y_staggered = base_label_y_upper_limit - (stagger_row_index - 1) * min_vertical_separation
    )
}


# --- Create the plot ---
plot <- ggplot() +
  # Tajima's D points
  geom_point(
    data = merged_data_filtered,
    aes(x = BIN_START, y = TajimaD, color = "Tajima's D"),
    size = 1.5, alpha = 0.8
  ) +
  # FST <= 0.1 points
  geom_point(
    data = merged_data_filtered %>% filter(FST_GROUP == "FST"),
    aes(x = BIN_START, y = FST_SCALED, color = "FST"),
    size = 1.5,
    alpha = 0.8
  ) +
  # FST > 0.1 points (highlighted in gold)
  geom_point(
    data = merged_data_filtered %>% filter(FST_GROUP == "FST > 0.1"),
    aes(x = BIN_START, y = FST_SCALED, color = "FST > 0.1"),
    size = 2.5, # Make them larger to stand out
    alpha = 1 # Make them fully opaque
  ) +
  # Centromere lines
  geom_vline(
    data = centromere_data_filtered,
    aes(xintercept = CENTROMERE_POS),
    linetype = "22",
    color = "gray30",
    size = 1.0,
    alpha=0.5,
  ) +
  # Add All Genes track (thinner, gray)
  geom_segment(data = genes_on_target_chr,
               aes(x = start, xend = end, y = gene_track_y_position_all_genes, yend = gene_track_y_position_all_genes_end, color = "All Genes"),
               size = 0.8, alpha = 0.6) # Thinner and slightly transparent


# Add Highlighted Gene Labels and Connecting Lines if there are genes to plot
if (nrow(genes_to_plot_final) > 0) {
  plot <- plot +
    # Add gene labels using ggrepel (only for highlighted genes)
    geom_text_repel(
      data = genes_to_plot_final,
      # Y-position for the segment line is 0 (the X-axis)
      aes(x = (start + end) / 2, y = -1, label = attribute),
      # Nudge label to its calculated staggered position
      nudge_y = genes_to_plot_final$label_y_staggered,
      color = "black",       # Text color
      segment.color = "black", segment.alpha = 0.25, # Connecting line color
      size = 6.0, # make gene name labels bigger 6-7 for chr10 5 for chr4
      fontface = "bold",
      box.padding = 0.4,          # Padding around label text
      point.padding = 0.4,        # Minimum distance from the data point
      max.overlaps = Inf,         # Ensures all labels are attempted to be plotted
      direction = "x",            # Prioritize horizontal repulsion once on assigned Y-level
      # Strictly confine labels to the Y = -5.0 to -2.0 range
      ylim = c(-5.0, -2.0),
      force = 100,                # Strong force to ensure horizontal separation
      min.segment.length = 0.0    # Ensure a line is always drawn
    )
} else {
  message(paste("No genes found for chromosome", target_chromosome, "that overlap highlighted FST windows."))
  message("This might be because your gene data is for a different chromosome or there are no significant FST peaks in the selected region for genes to overlap.")
}

# Add labels and scales
plot <- plot +
  scale_y_continuous(
    name = "Tajima's D",
    limits = c(-6.0, 5.0), # Overall y-axis limits (ensure labels are visible within this range)
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::label_number(),
    # Define secondary axis with custom transformation and breaks for 0-alignment
    sec.axis = sec_axis(~ . * 0.05, # New transformation: FST = TajimaD * 0.05 (0.05 is the slope)
                        name = "Weighted FST",
                        breaks = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.35), # New breaks to reflect the extended range
                        labels = scales::label_number()),
    expand = expansion(mult = c(0.05, 0.05)) # Still applies a small buffer *within* the set limits
  ) +
  scale_color_manual(
    values = c("Tajima's D" = "steelblue",
               "FST" = "#E69F00",
               "FST > 0.1" = "gold", # Set FST > 0.1 to gold
               "Genes" = "red", # Color for highlighted genes (used for connecting lines)
               "All Genes" = "orangered" # NEW: Color for the "All Genes" track
    ),
    name = ""
  )+
  scale_alpha_continuous(guide = "none") + # Hide alpha legend as alpha is set directly
  labs(
    title = paste0("Tajima's D & FST for Chr. ", target_chromosome, ": Testes Temp v Control"),
    x = "Genomic Position",
    color = "Metric"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 0.8))) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color = "steelblue", size = 20, face = "bold"),
    axis.text.y.right = element_text(color = "#E69F00", size = 20, face = "bold"),
    axis.text.x = element_text(color = "black", size = 20, face="bold"),
    axis.title.y = element_text(color = "steelblue", size = 20, face = "bold"),
    axis.title.y.right = element_text(color = "#E69F00", size = 20, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold", margin = margin(b = 5)),
    legend.text = element_text(face = "bold", size=20),
    legend.position = "top"
  )


# Save the plot as a high-resolution file
ggsave(paste0("TajimasD_and_FST_MT_MC_Chromosome_", target_chromosome, "_with_centromere_and_genes.png"), plot, width = 15, height = 9, dpi = 600)

# Print the plot
print(plot)
 

