# Title: GLM on Kimura data to analyse a shift in TE age based on treatment
# Author: Dr. Alice M. Godden

# Title: Plotting repeatlandscape from repeatmasker divsum files
# Author: Dr. Alice M. Godden

# step 1 - generate trimmed divsum files
# Load necessary libraries
library(dplyr)
library(readr)

# Define the directory containing the original divsum files
input_dir <- "~/Desktop/Projects/thermal stress/denovo_assembly_dnaseq_RM_Slow/TC"
output_dir <- "~/Desktop/Projects/thermal stress/denovo_assembly_dnaseq_RM_Slow/TC/cleaned_divsum"

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Get a list of all divsum files in the directory
divsum_files <- list.files(input_dir, pattern = "\\.divsum$", full.names = TRUE)


# Function to clean a divsum file
clean_divsum <- function(file_path, output_dir) {
  # Read the file
  lines <- readLines(file_path)
  
  # Find the line index where "Coverage for each repeat class and divergence (Kimura)" occurs
  start_index <- grep("Coverage for each repeat class and divergence \\(Kimura\\)", lines)
  
  # If the line is found, keep only the lines after it
  if (length(start_index) > 0) {
    cleaned_lines <- lines[(start_index + 1):length(lines)]
    
    # Write to a new file in the cleaned divsum directory
    output_file <- file.path(output_dir, basename(file_path))
    writeLines(cleaned_lines, output_file)
    
    print(paste("Cleaned file saved:", output_file))
  } else {
    print(paste("Skipping file (header not found):", file_path))
  }
}

# Apply the cleaning function to all divsum files
for (file in divsum_files) {
  clean_divsum(file, output_dir)
}

print("Cleaning process completed!")

# Load libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(ggsci)

# -----------------------------
# Define directories
dir_control <- "~/Desktop/Projects/thermal stress/denovo_assembly_dnaseq_RM_Slow/OC/cleaned_divsum"
dir_heat    <- "~/Desktop/Projects/thermal stress/denovo_assembly_dnaseq_RM_Slow/OT/cleaned_divsum"

# List all divsum files
files_control <- list.files(dir_control, pattern="\\.divsum$", full.names=TRUE)
files_heat    <- list.files(dir_heat, pattern="\\.divsum$", full.names=TRUE)

# Stop if no files found
if(length(files_control) == 0 | length(files_heat) == 0){
  stop("No .divsum files found in one of the directories!")
}

# -----------------------------
# Define genome size
genome_size <- 1373471384

# -----------------------------
# Function to process a file
process_divsum <- function(file, condition_label){
  df <- read_delim(file, delim = " ", col_types = cols())
  sample_name <- gsub(".*/|\\.divsum$", "", file)
  df$Sample <- sample_name
  df$Condition <- condition_label
  
  df_long <- pivot_longer(df, cols=-c(Div, Sample, Condition), names_to="Repeat_Type", values_to="Count")
  df_long <- df_long %>% mutate(Normalized = (Count / genome_size) * 100)
  return(df_long)
}

# -----------------------------
# Read all files and assign condition
df_control <- bind_rows(lapply(files_control, process_divsum, condition_label="Control"))
df_heat    <- bind_rows(lapply(files_heat, process_divsum, condition_label="Heat"))

# Combine
plot_df <- bind_rows(df_control, df_heat)

# -----------------------------
# Assign TE family groups
plot_df$Group <- ifelse(grepl("DNA", plot_df$Repeat_Type), "DNA",
                        ifelse(grepl("LINE", plot_df$Repeat_Type), "LINE",
                               ifelse(grepl("LTR", plot_df$Repeat_Type), "LTR",
                                      ifelse(grepl("SINE", plot_df$Repeat_Type), "SINE",
                                             ifelse(grepl("RC", plot_df$Repeat_Type), "RC",
                                                    ifelse(grepl("Satellite", plot_df$Repeat_Type), "Satellite",
                                                           ifelse(grepl("Retrotransposon", plot_df$Repeat_Type), "Retrotransposon", NA)))))))

# Remove NA groups
plot_df <- plot_df %>% filter(!is.na(Group))

# -----------------------------
# Plot: facet by TE family, color by condition
ggplot(plot_df, aes(x = Div, y = Normalized, color = Condition, group = interaction(Sample, Repeat_Type))) +
  geom_line(size = 0.8, alpha = 0.4) +
  scale_color_manual(values = c("Control" = "blue", "Heat" = "orange")) +
  coord_cartesian(xlim = c(0,30), ylim = c(0,0.8)) +
  facet_wrap(~ Group, scales = "free_y") +
  theme_classic() +
  xlab("Kimura Substitution Level (CpG adjusted)") +
  ylab("Percent of Genome (%)") +
  labs(color = "Condition") +
  ggtitle("Repeat Landscapes by TE Family and Condition (Testes)") +
  theme(
    axis.text = element_text(size = 11, face = "bold", color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.text = element_text(face = "bold")
  )

# free y axis
# Plot: facet by TE family, color by condition
ggplot(plot_df, aes(x = Div, y = Normalized, color = Condition, group = interaction(Sample, Repeat_Type))) +
  geom_line(size = 0.8, alpha = 0.4) +
  scale_color_manual(values = c("Control" = "blue", "Heat" = "orange")) +
  coord_cartesian(xlim = c(0,30)) +   # Remove ylim to auto-scale
  facet_wrap(~ Group, scales = "free_y") +
  theme_classic() +
  xlab("Kimura Substitution Level (CpG adjusted)") +
  ylab("Percent of Genome (%)") +
  labs(color = "Condition") +
  ggtitle("Repeat Landscapes by TE Family and Condition (Ovaries)") +
  theme(
    axis.text = element_text(size = 11, face = "bold", color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.text = element_text(face = "bold")
  )

ggsave("Repeat_Landscapes_Ovaries.png", width = 8, height = 8, dpi = 300)


# stats
# ==============================================================
# Repeat Landscape & TE Shift Analysis (testes)
# ==============================================================

# Load libraries
library(tidyverse)
library(MASS)  # Negative Binomial
library(pscl)  # Zero-inflated models

# -----------------------------
# Step 1: Define directories
dir_control <- "~/Desktop/Projects/thermal stress/denovo_assembly_dnaseq_RM_Slow/OC/cleaned_divsum"
dir_heat    <- "~/Desktop/Projects/thermal stress/denovo_assembly_dnaseq_RM_Slow/OT/cleaned_divsum"

# List files
files_control <- list.files(dir_control, pattern="\\.divsum$", full.names=TRUE)
files_heat    <- list.files(dir_heat, pattern="\\.divsum$", full.names=TRUE)

if(length(files_control) == 0 | length(files_heat) == 0){
  stop("No .divsum files found in one of the directories!")
}

# -----------------------------
# Step 2: Define genome size
genome_size <- 1373471384

# -----------------------------
# Step 3: Function to read and process divsum files
process_divsum <- function(file, condition_label){
  df <- read_delim(file, delim = " ", col_types = cols())
  sample_name <- gsub(".*/|\\.divsum$", "", file)
  df$Sample <- sample_name
  df$Condition <- condition_label
  
  df_long <- pivot_longer(df, cols = -c(Div, Sample, Condition), 
                          names_to = "Repeat_Type", values_to = "Count")
  df_long <- df_long %>% mutate(Count = ifelse(Count < 0, 0, Count),
                                Normalized = (Count / genome_size) * 100)
  return(df_long)
}

# -----------------------------
# Step 4: Read all files
df_control <- bind_rows(lapply(files_control, process_divsum, condition_label="Control"))
df_heat    <- bind_rows(lapply(files_heat, process_divsum, condition_label="Heat"))

plot_df <- bind_rows(df_control, df_heat)

# -----------------------------
# Step 5: Assign TE Families
plot_df <- plot_df %>%
  mutate(TE_Family = case_when(
    str_detect(Repeat_Type, "DNA")        ~ "DNA",
    str_detect(Repeat_Type, "LINE")       ~ "LINE",
    str_detect(Repeat_Type, "LTR")        ~ "LTR",
    str_detect(Repeat_Type, "SINE")       ~ "SINE",
    str_detect(Repeat_Type, "RC")         ~ "RC",
    str_detect(Repeat_Type, "Satellite")  ~ "Satellite",
    str_detect(Repeat_Type, "Retrotransposon") ~ "Retrotransposon",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(TE_Family))

plot_df$Condition <- as.factor(plot_df$Condition)
plot_df$TE_Family <- as.factor(plot_df$TE_Family)

# -----------------------------
# Step 6: Overdispersion function
overdispersion_test <- function(model) {
  rdf <- df.residual(model)
  res_dev <- sum(residuals(model, type="pearson")^2)
  return(res_dev / rdf)
}

# -----------------------------
# Step 7: Run GLM / NB / ZIP / ZINB models per TE family
shift_results <- tibble()

for(te in unique(plot_df$TE_Family)){
  cat("\n📌 Testing TE Family:", te, "...\n")
  
  te_data <- filter(plot_df, TE_Family == te)
  
  # Fit models
  poisson_model <- tryCatch(glm(Count ~ Condition + Div, data=te_data, family=poisson), error=function(e) NULL)
  nb_model      <- tryCatch(glm.nb(Count ~ Condition + Div, data=te_data, control=glm.control(maxit=25)), error=function(e) NULL)
  zip_model     <- tryCatch(zeroinfl(Count ~ Condition + Div | 1, data=te_data, dist="poisson"), error=function(e) NULL)
  zinb_model    <- tryCatch(zeroinfl(Count ~ Condition + Div | 1, data=te_data, dist="negbin"), error=function(e) NULL)
  
  models <- list(
    Poisson = poisson_model,
    NegBin  = nb_model,
    ZIP     = zip_model,
    ZINB    = zinb_model
  )
  
  for(mname in names(models)){
    mod <- models[[mname]]
    if(!is.null(mod)){
      coef_sum <- summary(mod)$coefficients
      pval <- if("ConditionHeat" %in% rownames(coef_sum)) coef_sum["ConditionHeat","Pr(>|z|)"] else NA
      est  <- if("ConditionHeat" %in% rownames(coef_sum)) coef_sum["ConditionHeat","Estimate"] else NA
      shift <- if(!is.na(pval) & pval<0.05 & est>0) "Increase" else if(!is.na(pval) & pval<0.05 & est<0) "Decrease" else "No Significant Shift"
      star <- if(!is.na(pval) & pval<0.001) "***" else if(!is.na(pval) & pval<0.01) "**" else if(!is.na(pval) & pval<0.05) "*" else "ns"
      overdisp <- if(mname %in% c("NegBin","ZINB")) overdispersion_test(mod) else NA
      
      shift_results <- bind_rows(shift_results, tibble(
        TE_Family = te,
        Model = mname,
        AIC = AIC(mod),
        Estimate = est,
        P_Value = pval,
        Overdispersion = overdisp,
        Shift = shift,
        Significance = star
      ))
    }
  }
}

# -----------------------------
# Step 8: Save results
write.csv(shift_results, file.path("~/Desktop/Projects/thermal stress/denovo_assembly_dnaseq_RM_Slow", 
                                   "TE_Family_Shift_Results_Ovaries.csv"), row.names=FALSE)
print(shift_results)

cat("\n✅ Analysis complete! Results saved.\n")


