# Load necessary libraries
library(ggplot2)
library(viridis)
library(gridExtra)
library(dplyr)

# Define the genome sizes for each dataset
genomes_sizes <- c(39199570013, 40253710347, 41736660246, 44745217617)
titles <- c("Ovaries Control", "Ovaries Temperature", "Testes Control", "Testes Temperature")

# Function to process each input file and create a pie chart with percentage labels
process_input_file <- function(file_path, title, genomes_size) {
  # Read the input file
  KimuraDistance <- read.csv(file_path, sep = " ")
  
  # Calculate normalized values
  kd_melt <- melt(KimuraDistance, id = "Div")
  kd_melt$norm <- kd_melt$value / genomes_size * 100
  
  # Create a new column based on substring matches
  kd_melt$group <- ifelse(grepl("DNA", kd_melt$variable), "DNA",
                          ifelse(grepl("LINE", kd_melt$variable), "LINE",
                                 ifelse(grepl("LTR", kd_melt$variable), "LTR",
                                        ifelse(grepl("SINE", kd_melt$variable), "SINE",
                                               ifelse(grepl("RC", kd_melt$variable), "RC",
                                                      ifelse(grepl("Satellite", kd_melt$variable), "Satellite",
                                                             ifelse(grepl("Retrotransposon", kd_melt$variable), "Retrotransposon", NA)))))))
  
  # Filter out NA rows
  kd_melt_filtered <- kd_melt[!is.na(kd_melt$group), ]
  
  # Summarize data for pie chart
  pie_data <- aggregate(norm ~ group, data = kd_melt_filtered, sum)
  
  # Calculate percentage labels
  pie_data$label <- paste0(pie_data$group, "\n", round(pie_data$norm, 1), "%")
  
  # Calculate "Unmasked" category to make pie chart add up to 100%
  total <- sum(pie_data$norm)
  pie_data <- rbind(pie_data, c("Unmasked", 100 - total))
  
  # Calculate positions for each segment
  pie_data <- pie_data %>%
    arrange(desc(norm)) %>%
    mutate(end = cumsum(norm), start = lag(end, default = 0))
  
  # Create custom fill colors
  custom_colors <- c(DNA = "#000004FF", LINE = "#3B0F70FF", LTR = "#8C2981FF",
                     SINE = "#FE9F6DFF", RC = "#DE4968FF", Satellite ="#FCFDBFFF",
                     Retrotransposon = "white", Unmasked = "grey50")
  
  
  # Create pie chart with percentage labels and custom color scheme
  pie_chart <- ggplot(pie_data, aes(x = 1, y = norm, fill = group)) +
    geom_rect(aes(xmin = -1, xmax = 1, ymin = c(0, head(end, -1)), ymax = end, fill = group),
              color =FALSE) +
    geom_text(aes(label = label, y = (start + end) / 2),
              color = "darkgrey",
              size = 3.,
              fontface = "bold") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = custom_colors) +  # Use custom fill scale
    labs(fill = "") +
    theme_void() +
    ggtitle(paste("RNA-seq:", title, "Genome Composition"))
  
  return(pie_chart)
}

# Create a list to store individual pie chart objects
file_paths <- c("FC_all.divsum.csv", "FT_all.divsum.csv", "MC_all.divsum.csv", "MT_all.divsum.csv")
plots <- lapply(seq_along(file_paths), function(i) {
  process_input_file(file_paths[i], titles[i], genomes_sizes[i])
})

# Combine plots using grid.arrange
final_plot <- grid.arrange(grobs = plots, ncol = 2)  # Adjust ncol as needed

# Show the final plot
final_plot



----
# to get % composition of TES
  
# Load necessary libraries
library(ggplot2)
library(viridis)
library(gridExtra)

# Define the genome sizes for each dataset
genomes_sizes <- c(39199570013, 40253710347, 41736660246, 44745217617)
titles <- c("Ovaries Control", "Ovaries Temperature", "Testes Control", "Testes Temperature")

# Function to process each input file and create a pie chart with percentage labels
process_input_file <- function(file_path, title, genomes_size) {
  # Read the input file
  KimuraDistance <- read.csv(file_path, sep = " ")
  
  # Calculate normalized values
  kd_melt <- melt(KimuraDistance, id = "Div")
  kd_melt$norm <- kd_melt$value / genomes_size * 100
  
  # Create a new column based on substring matches
  kd_melt$group <- ifelse(grepl("DNA", kd_melt$variable), "DNA",
                          ifelse(grepl("LINE", kd_melt$variable), "LINE",
                                 ifelse(grepl("LTR", kd_melt$variable), "LTR",
                                        ifelse(grepl("SINE", kd_melt$variable), "SINE",
                                               ifelse(grepl("RC", kd_melt$variable), "RC",
                                                      ifelse(grepl("Satellite", kd_melt$variable), "Satellite",
                                                             ifelse(grepl("Retrotransposon", kd_melt$variable), "Retrotransposon", NA)))))))
  
  # Filter out NA rows
  kd_melt_filtered <- kd_melt[!is.na(kd_melt$group), ]
  
  # Summarize data for pie chart
  pie_data <- aggregate(norm ~ group, data = kd_melt_filtered, sum)
  
  # Calculate percentage labels
  pie_data$label <- paste0(pie_data$group, "\n", round(pie_data$norm, 1), "%")
  
  # Calculate positions for each segment
  pie_data <- pie_data %>%
    arrange(desc(norm)) %>%
    mutate(end = cumsum(norm), start = lag(end, default = 0))
  
  # Create pie chart with percentage labels and magma color scheme
  pie_chart <- ggplot(pie_data, aes(x = 1, y = norm, fill = group)) +
    geom_rect(aes(x = 1, xmin = 0, xmax = 1, ymin = start, ymax = end), color = "white") +
    geom_text(aes(x = 1, y = (start + end) / 2, label = label),
              color="darkgrey",
              size=3.,
              fontface = "bold") +
    coord_polar(theta = "y", start = 0) +
    scale_fill_viridis_d(option = "magma") +  # Use magma color scheme
    labs(fill = "") +
    theme_void() +
    ggtitle(paste("RNA-seq:", title, "Genome Composition"))
  
  return(pie_chart)
}

# Create a list to store individual pie chart objects
plots <- lapply(seq_along(file_paths), function(i) {
  process_input_file(file_paths[i], titles[i], genomes_sizes[i])
})

# Combine plots using grid.arrange
final_plot <- grid.arrange(grobs = plots, ncol = 2)  # Adjust ncol as needed

# Show the final plot
final_plot


#### without percentages

# Load necessary libraries
library(ggplot2)
library(viridis)
library(gridExtra)
library(dplyr)

# Define the genome sizes for each dataset
genomes_sizes <- c(39199570013, 40253710347, 41736660246, 44745217617)
titles <- c("Ovaries Control", "Ovaries Temperature", "Testes Control", "Testes Temperature")

# Function to process each input file and create a pie chart with percentage labels
process_input_file <- function(file_path, title, genomes_size) {
  # Read the input file
  KimuraDistance <- read.csv(file_path, sep = " ")
  
  # Calculate normalized values
  kd_melt <- melt(KimuraDistance, id = "Div")
  kd_melt$norm <- kd_melt$value / genomes_size * 100
  
  # Create a new column based on substring matches
  kd_melt$group <- ifelse(grepl("DNA", kd_melt$variable), "DNA",
                          ifelse(grepl("LINE", kd_melt$variable), "LINE",
                                 ifelse(grepl("LTR", kd_melt$variable), "LTR",
                                        ifelse(grepl("SINE", kd_melt$variable), "SINE",
                                               ifelse(grepl("RC", kd_melt$variable), "RC",
                                                      ifelse(grepl("Satellite", kd_melt$variable), "Satellite",
                                                             ifelse(grepl("Retrotransposon", kd_melt$variable), "Retrotransposon", NA)))))))
  
  # Filter out NA rows
  kd_melt_filtered <- kd_melt[!is.na(kd_melt$group), ]
  
  # Summarize data for pie chart
  pie_data <- aggregate(norm ~ group, data = kd_melt_filtered, sum)
  
  # Calculate "Unmasked" category to make pie chart add up to 100%
  total <- sum(pie_data$norm)
  pie_data <- rbind(pie_data, c("Unmasked", 100 - total))
  
  # Calculate positions for each segment
  pie_data <- pie_data %>%
    arrange(desc(norm)) %>%
    mutate(end = cumsum(norm), start = lag(end, default = 0))
  
  # Create custom fill colors
  custom_colors <- c(DNA = "#000004FF", LINE = "#3B0F70FF", LTR = "#8C2981FF",
                     SINE = "#FE9F6DFF", RC = "#DE4968FF", Satellite ="#FCFDBFFF",
                     Retrotransposon = "white", Unmasked = "grey50")
  
  # Create pie chart without percentage labels
  pie_chart <- ggplot(pie_data, aes(x = 1, y = norm, fill = group)) +
    geom_rect(aes(xmin = -1, xmax = 1, ymin = c(0, head(end, -1)), ymax = end, fill = group),
              color = FALSE) +
    geom_text(aes(label = ""),  # Empty string to remove labels
              color = "darkgrey",
              size = 3.,
              fontface = "bold") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = custom_colors) +  # Use custom fill scale
    labs(fill = "") +  # Remove legend labels
    theme_void() +
    theme(
      legend.text = element_text(face = "bold"),  # Make legend text bold
      plot.title = element_text(face = "bold")  # Make plot title bold
    ) +
    ggtitle(paste("RNA-seq:", title, "Genome Composition"))
  
  return(pie_chart)
}

# Create a list to store individual pie chart objects
file_paths <- c("FC_all.divsum.csv", "FT_all.divsum.csv", "MC_all.divsum.csv", "MT_all.divsum.csv")
plots <- lapply(seq_along(file_paths), function(i) {
  process_input_file(file_paths[i], titles[i], genomes_sizes[i])
})

# Combine plots using grid.arrange
final_plot <- grid.arrange(grobs = plots, ncol = 2)  # Adjust ncol as needed

# Show the final plot
final_plot

