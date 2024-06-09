#Title: Telescope output read in to make final counts file
#Author: Dr. Alice M. Godden

# Load necessary libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

# Function to read in a TSV file and add the file name as a column
customized_read_tsv <- function(file) {
  read_tsv(file, col_types = cols(
    transcript = col_character(),
    count = col_double()
  )) %>%
    mutate(sample = gsub(".*telescope-TE_counts_(SAMN[0-9]+)\\.tsv", "\\1", file))
}

# Directory containing the TSV files
file_dir <- "~/Desktop/projectbear/newgenome/rnaseq/Telescope/telescope_output"

# List all the TSV files in the directory
f_files <- list.files(file_dir, pattern = "*.tsv", full.names = TRUE)

# Read and process each TSV file
df_list <- lapply(f_files, customized_read_tsv)

# Combine all data frames into a single data frame
merged <- bind_rows(df_list)

# Pivot wider to have each sample as a column
final_merged <- merged %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = list(count = 0))

# Write the merged data to a TSV file
write_tsv(final_merged, file = "telescope_merged_counts.tsv")

# Print the first few rows of the merged data
print(head(final_merged))
