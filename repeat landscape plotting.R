#############################
###Plot Kimura distance######
#############################

library(reshape)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(gridExtra)

sessionInfo()
#R version 3.6.1 (2019-07-05)
#Platform: x86_64-conda_cos6-linux-gnu (64-bit)
#Running under: Ubuntu 20.04.1 LTS

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
#other attached packages:
#  [1] gridExtra_2.3     forcats_0.5.0     stringr_1.4.0     dplyr_0.8.5       purrr_0.3.4      
#  [6] readr_1.3.1       tidyr_1.1.0       tibble_3.0.1      tidyverse_1.3.0   hrbrthemes_0.8.0 
#  [11] viridis_0.5.1     viridisLite_0.3.0 ggplot2_3.3.0     reshape_0.8.8    

KimuraDistance <- read.csv("FC_all.divsum.csv",sep=" ")

#add here the genome size in bp
genomes_size=39199570013
  # ensembl 1x genome 1373471384
  # FC merged = 39199570013
  # FT merged = 40253710347
  # MC merged = 41736660246
  # MT merged = 44745217617 

kd_melt = melt(KimuraDistance,id="Div")
kd_melt$norm = kd_melt$value/genomes_size * 100




# Create a new column in kd_melt based on the substring match
kd_melt$group <- ifelse(grepl("DNA", kd_melt$variable), "DNA",
                        ifelse(grepl("LINE", kd_melt$variable), "LINE",
                               ifelse(grepl("LTR", kd_melt$variable), "LTR",
                                      ifelse(grepl("SINE", kd_melt$variable), "SINE",
                                             ifelse(grepl("RC", kd_melt$variable), "RC",
                                                    ifelse(grepl("Satellite", kd_melt$variable), "Satellite",
                                                           ifelse(grepl("Retrotransposon", kd_melt$variable), "Retrotransposon", NA)))))))

# Define a custom color palette using Viridis "magma" color scheme
group_colors <- magma(length(unique(kd_melt$group)))

# Assign colors to each group based on the custom palette
names(group_colors) <- unique(kd_melt$group)



# to remove NA

library(viridis)

# Filter out rows where group is NA
kd_melt_filtered <- kd_melt[!is.na(kd_melt$group), ]

# Define a custom color palette using Viridis "magma" color scheme
group_colors <- magma(length(unique(kd_melt_filtered$group)))

# Assign colors to each group based on the custom palette
names(group_colors) <- unique(kd_melt_filtered$group)

# Use the custom color palette in ggplot
ggplot(kd_melt_filtered, aes(fill = group, y = norm, x = Div)) + 
  geom_bar(position = "stack", stat = "identity", color = FALSE) +
  scale_fill_manual(values = group_colors) +  # Use manual fill colors
  theme_classic() +
  xlab("Kimura substitution level") +
  ylab("Percent of the genome") + 
  labs(fill = "") +
  coord_cartesian(xlim = c(0, 55)) +
  ggtitle("RNA-seq Ovaries control: Kimura Substitution Level") +
  theme(axis.text = element_text(size = 11, face = "bold", color = "black"), 
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold"))

# how to have all four plots together in a facet grid

file_paths <- c("FC_all.divsum.csv", "FT_all.divsum.csv", "MC_all.divsum.csv", "MT_all.divsum.csv")
titles <- c("Ovaries Control", "Ovaries Temperature", "Testes Control", "Testes Temperature")
genomes_sizes <- c(39199570013, 40253710347, 41736660246, 44745217617)


library(ggplot2)
library(viridis)

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
  
  # Define custom color palette
  group_colors <- magma(length(unique(kd_melt_filtered$group)))
  names(group_colors) <- unique(kd_melt_filtered$group)
  
  # Create ggplot object
  gg <- ggplot(kd_melt_filtered, aes(fill = group, y = norm, x = Div)) + 
    geom_bar(position = "stack", stat = "identity", color = FALSE) +
    scale_fill_manual(values = group_colors) +
    theme_classic() +
    xlab("Kimura substitution level") +
    ylab("Percent of the genome") + 
    labs(fill = "") +
    coord_cartesian(xlim = c(0, 55)) +
    ggtitle(paste("RNA-seq:", title, "Kimura Substitution Level")) +
    theme(axis.text = element_text(size = 11, face = "bold", color = "black"), 
          axis.title = element_text(size = 12, face = "bold"),
          # Add legend.text to set legend text properties
          legend.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))
  
  return(gg)
}

library(gridExtra)

# Create a list to store individual ggplot objects
plots <- lapply(seq_along(file_paths), function(i) {
  process_input_file(file_paths[i], titles[i], genomes_sizes[i])
})

# Combine plots using facet_wrap
final_plot <- do.call(gridExtra::grid.arrange, c(plots, ncol = 2))  # Adjust ncol as needed

# Show the final plot
final_plot

#To control the y axes limits for each plot
# Load necessary libraries
library(reshape)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(gridExtra)

# Print session information for debugging purposes
sessionInfo()

# Define a function to process input files and generate plots
process_input_file <- function(file_path, title, genomes_size, y_limit) {
  # Read the input file
  KimuraDistance <- read.csv(file_path, sep = " ")
  
  # Melt the data frame
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
  
  # Define a custom color palette using Viridis "magma" color scheme
  group_colors <- magma(length(unique(kd_melt_filtered$group)))
  names(group_colors) <- unique(kd_melt_filtered$group)
  
  # Create ggplot object
  gg <- ggplot(kd_melt_filtered, aes(fill = group, y = norm, x = Div)) + 
    geom_bar(position = "stack", stat = "identity", color = FALSE) +
    scale_fill_manual(values = group_colors) +
    theme_classic() +
    xlab("Kimura substitution level (CpG adjusted)") +
    ylab("Percent of the genome (%)") + 
    labs(fill = "") +
    coord_cartesian(xlim = c(0, 30), ylim = c(0, y_limit)) +
    ggtitle(paste("RNA-seq:", title, "Kimura Substitution Level")) +
    theme(
      axis.text = element_text(size = 11, face = "bold", color = "black"), 
      axis.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      panel.grid.major = element_line(size = 0.1, linetype = 'solid', color = "gray"),  # Customize major gridlines
      panel.grid.minor = element_line(size = 0.05, linetype = 'solid', color = "gray")  # Customize minor gridlines
    )
  
  return(gg)
}

# File paths, titles, genome sizes, and y-axis limits
file_paths <- c("FC_all.divsum.csv", "FT_all.divsum.csv", "MC_all.divsum.csv", "MT_all.divsum.csv")
titles <- c("Ovaries Control", "Ovaries Temperature", "Testes Control", "Testes Temperature")
genomes_sizes <- c(39199570013, 40253710347, 41736660246, 44745217617)
y_limits <- c(1.7, 1.7, 7, 7)  # y-axis limits for each plot

# Create a list to store individual ggplot objects
plots <- lapply(seq_along(file_paths), function(i) {
  process_input_file(file_paths[i], titles[i], genomes_sizes[i], y_limits[i])
})

# Combine plots using grid.arrange
final_plot <- do.call(gridExtra::grid.arrange, c(plots, ncol = 2))  # Adjust ncol as needed

# Show the final plot
final_plot


