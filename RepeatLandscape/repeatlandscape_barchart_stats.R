#Title: Plotting bar charts with stats for repeat landscape divsum files
#Author : Dr. Alice M. Godden

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(readr)
library(viridis)

# Load the control and treatment divsum files
control <- read_delim("MC_all.divsum.csv", delim = " ", col_types = cols())
treatment <- read_delim("MT_all.divsum.csv", delim = " ", col_types = cols())

# Add group labels
control <- control %>% mutate(Group = "Control")
treatment <- treatment %>% mutate(Group = "Temperature")  # Changed Treatment to Temperature

# Combine datasets
data <- bind_rows(control, treatment)

# Check column names to confirm TE categories
print(colnames(data))

# Dynamically select and group relevant TE families
data_long <- data %>%
  select(Group, contains("DNA"), contains("LINE"), contains("LTR"), contains("RC"), contains("Satellite"), contains("SINE")) %>%
  pivot_longer(cols = -Group, names_to = "Subcategory", values_to = "Count") %>%
  mutate(Family = case_when(
    grepl("DNA", Subcategory) ~ "DNA",
    grepl("LINE", Subcategory) ~ "LINE",
    grepl("LTR", Subcategory) ~ "LTR",
    grepl("RC", Subcategory) ~ "RC",
    grepl("Satellite", Subcategory) ~ "Satellite",
    grepl("SINE", Subcategory) ~ "SINE",
    TRUE ~ "Other"
  ))

# Summarize data by family and group
family_data <- data_long %>%
  group_by(Family, Group) %>%
  summarise(Total_Count = sum(Count, na.rm = TRUE), .groups = "drop")

# Perform statistical tests for each family
stat_results <- data_long %>%
  group_by(Family) %>%
  summarise(
    is_normal = ifelse(
      shapiro.test(Count[Group == "Control"])$p.value > 0.05 &
        shapiro.test(Count[Group == "Temperature"])$p.value > 0.05, TRUE, FALSE
    ),
    p_value = ifelse(
      is_normal,
      t.test(Count ~ Group)$p.value,
      wilcox.test(Count ~ Group)$p.value
    ),
    test_used = ifelse(
      is_normal,
      "T-test (Normality confirmed for both groups)",
      "Wilcoxon test (Non-normal data)"
    )
  )

# Print which statistical test was used and why
print("Statistical tests used:")
print(stat_results %>% select(Family, test_used, p_value))

# Add significance levels based on p-values
stat_results <- stat_results %>%
  mutate(
    significance = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01  ~ "**",
      p_value <= 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Merge statistical results back to the summarized family data
family_data <- left_join(family_data, stat_results, by = "Family")

# Bar chart with significance stars
ggplot(family_data, aes(x = Family, y = Total_Count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", fontface = "bold") +
  geom_text(data = stat_results, aes(x = Family, y = max(family_data$Total_Count) * 1.1, label = significance), inherit.aes = FALSE) +
  labs(
    title = "Testes Repeat Landscape",
    x = "Family",
    y = "Total Counts",
    fill = "Group"
  ) +
  scale_fill_manual(
    values = c("Control" = "#5C1D6E", "Temperature" = "#E84C33")  # Specify custom colors
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face="bold")) +
  theme(
    plot.title = element_text(face = "bold"),  # Title
    axis.title = element_text(face = "bold", color ="black"),  # Axis titles
    axis.text = element_text(face = "bold", color = "black"),   # Axis text
    legend.title = element_text(face = "bold", color = "black"),  # Legend title
    legend.text = element_text(face = "bold", color = "black")   # Legend text
  )

# Save the plot
ggsave("TESTES_te_family_bar_chart_with_temperature.png", width = 10, height = 6)
