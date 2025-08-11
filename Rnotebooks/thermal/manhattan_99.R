# Title : FST manhattan 99th percentile
# Author : Dr. Alice M. Godden

# --- IMPORTANT: Run all library() calls first! ---
  # Load necessary packages
  library(dplyr)   # For data manipulation (e.g., filter, %>%, mutate, group_by, summarise)
library(ggplot2) # For creating the Manhattan plot

# --- Load your data ---
# This reads your FST data into the 'fst3' data frame.
fst3 <- read.table("fst_outputMTvMC.windowed.weir.fst", sep = '\t', header = TRUE)

# --- Data Preparation for fst3 ---
cat("--- fst3 Data Preparation ---\n")
cat("Initial dimensions of fst3:", dim(fst3)[1], "rows, ", dim(fst3)[2], "columns\n")
cat("Initial unique values in CHROM column:", unique(fst3$CHROM), "\n")


# 1. Filter out rows where CHROM is "MT" and convert CHROM to numeric
fst3 <- fst3 %>%
  filter(CHROM != "MT") %>%
  mutate(CHROM = as.numeric(as.character(CHROM))) # Convert to character first to handle non-numeric if any, then to numeric

cat("Dimensions after filtering 'MT' chromosome:", dim(fst3)[1], "rows\n")
cat("Unique values in CHROM column after filtering and numeric conversion:", unique(fst3$CHROM), "\n")


# 2. Create SNP column (if not already present and needed for identification)
length_ <- dim(fst3)[1] # Get the updated number of rows after filtering
fst3$SNP <- paste('SNP', 1:length_)

cat("First few rows of fst3 after all modifications:\n")
print(head(fst3))
cat("\n")

# 99th percentile calculation for highlighting
fst_99th_percentile <- quantile(fst3$WEIGHTED_FST, 0.99, na.rm = TRUE)
cat("99th percentile of WEIGHTED_FST:", fst_99th_percentile, "\n")

# Identify the SNPs that have a WEIGHTED_FST value above the 99th percentile.
# We'll use this to create a 'highlight' column for ggplot2.
fst3 <- fst3 %>%
  mutate(is_highlight = ifelse(WEIGHTED_FST > fst_99th_percentile, TRUE, FALSE))

cat("Number of SNPs to be highlighted (above 99th percentile):", sum(fst3$is_highlight), "\n\n")


# --- Prepare data for ggplot2 Manhattan plot ---
# Calculate cumulative positions for the x-axis
data_plot <- fst3 %>%
  # Group by chromosome and calculate the maximum base pair position for each chromosome
  group_by(CHROM) %>%
  summarise(chr_len = max(BIN_START)) %>%
  # Calculate cumulative length for each chromosome
  mutate(total_len = cumsum(as.numeric(chr_len)) - chr_len) %>%
  # Join back to the original data to get the adjusted position for each SNP
  select(CHROM, total_len) %>%
  left_join(fst3, ., by = c("CHROM")) %>%
  mutate(BP_adjusted = BIN_START + total_len) # Calculate the adjusted position for plotting


# Prepare chromosome midpoints for x-axis labels
axis_df <- data_plot %>%
  group_by(CHROM) %>%
  summarise(center = mean(BP_adjusted))


# --- Manhattan Plot Creation using ggplot2 ---
manhattan_plot <- ggplot(data_plot, aes(x = BP_adjusted, y = WEIGHTED_FST)) +
  # Add points for all SNPs, colored by chromosome
  geom_point(aes(color = as.factor(CHROM)), alpha = 0.8, size = 1.3) +
  # Use specific colors for alternating chromosomes
  scale_color_manual(values = rep(c("lightblue", "navy"), unique(length(axis_df$CHROM)))) +
  # Add a separate layer for highlighted SNPs in gold
  geom_point(data = filter(data_plot, is_highlight == TRUE),
             color = "goldenrod",
             size = 2.5,
             pch = 20, # Use filled circle for highlighted points
             alpha = 1) +
  # Set y-axis limits
  ylim(0.0, 0.7) +
  # Add labels and title
  labs(y = 'WEIGHTED FST',
       x = 'CHROMOSOME',
       title = 'FST Testes Temp v Control windowed 100kb') +
  # Set x-axis breaks and labels to show chromosome numbers
  scale_x_continuous(label = axis_df$CHROM, breaks = axis_df$center) +
  # Apply the desired theme elements (bold text)
  theme_minimal() + # Use a minimal theme as a base
  theme(
    legend.position = "none", # Hide the chromosome color legend
    panel.grid.major.x = element_blank(), # Remove vertical grid lines
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(face = "bold", size = 14, colour = "black"), # Bold x-axis tick labels
    axis.text.y = element_text(face = "bold", size = 14, colour = "black"), # Bold y-axis tick labels
    axis.title.x = element_text(face = "bold", size = 14, colour = "black", margin = margin(t = 10)), # Bold x-axis title
    axis.title.y = element_text(face = "bold", size = 14, colour = "black", margin = margin(r = 10)), # Bold y-axis title
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18) # Bold and center plot title
  )

# Print the plot to the plots window
print(manhattan_plot) # this don't work- don't know why. 

# --- Save the plot ---
# Save the ggplot object directly using ggsave
ggsave('manhattan_plot_MTvMC_HS.png', plot = manhattan_plot,
       width = 9, height = 5, units = "in", dpi = 600) # Save with higher resolution

cat("\nManhattan plot saved as 'manhattan_plot_Fposter_F3SvF3L.png' in your working directory.\n")
