# Title: GLM on Kimura data to analyse a shift in TE age based on treatment
# Author: Dr. Alice M. Godden

# GLM selection

# Load necessary libraries
library(tidyverse)
library(MASS)  # For Negative Binomial Model (glm.nb)
library(pscl)  # For Zero-Inflated Model (zeroinfl)

# Step 1: Read Data
control_data <- read.delim("FC_all.divsum.csv", sep="", header=FALSE, check.names=FALSE, stringsAsFactors=FALSE)
temperature_data <- read.delim("FT_all.divsum.csv", sep="", header=FALSE, check.names=FALSE, stringsAsFactors=FALSE)

# Step 2: Assign Column Names
colnames(control_data) <- as.character(control_data[1, ])
colnames(temperature_data) <- as.character(temperature_data[1, ])
control_data <- control_data[-1, ]
temperature_data <- temperature_data[-1, ]
rownames(control_data) <- NULL
rownames(temperature_data) <- NULL

# Step 3: Convert Columns to Numeric & Remove Negatives
convert_numeric <- function(df) {
  df %>% mutate(across(where(is.character), ~ suppressWarnings(as.numeric(.)))) %>%
    replace(is.na(.), 0)  # Replace NA with 0
}

control_data <- convert_numeric(control_data)
temperature_data <- convert_numeric(temperature_data)
control_data[control_data < 0] <- 0
temperature_data[temperature_data < 0] <- 0
cat("✅ No negative values remain.\n")

# Step 4: Add Condition Column & Merge Data
control_data$Condition <- "Control"
temperature_data$Condition <- "Heat_Stress"
merged_data <- bind_rows(control_data, temperature_data)

# Step 5: Reshape Data to Long Format
data_long <- merged_data %>%
  pivot_longer(cols = -c(Div, Condition), names_to = "Repeat_Family", values_to = "Count")

# Step 6: Assign TE Families
data_long <- data_long %>%
  mutate(TE_Family = case_when(
    str_detect(Repeat_Family, "DNA")  ~ "DNA",
    str_detect(Repeat_Family, "LINE") ~ "LINE",
    str_detect(Repeat_Family, "LTR")  ~ "LTR",
    str_detect(Repeat_Family, "SINE") ~ "SINE",
    str_detect(Repeat_Family, "RC")   ~ "RC",
    str_detect(Repeat_Family, "Satellite")  ~ "Satellite",
    TRUE ~ "Other"
  ))

data_long$Condition <- as.factor(data_long$Condition)
data_long$TE_Family <- as.factor(data_long$TE_Family)

# Step 7: Remove Any NA Values in TE_Family
data_long <- data_long %>% filter(!is.na(TE_Family))

# Step 8: Run Separate Models for Each TE Family and Test for Shift
shift_results <- tibble()

# Overdispersion Function
overdispersion_test <- function(model) {
  rdf <- df.residual(model)  # Residual degrees of freedom
  res_dev <- sum(residuals(model, type = "pearson")^2)  # Pearson Residuals
  return(res_dev / rdf)  # Overdispersion parameter
}

for (te in unique(data_long$TE_Family)) {
  
  cat("\n📌 Testing TE Family:", te, "...\n")
  
  te_data <- data_long %>% filter(TE_Family == te)
  
  # Fit Poisson Model
  poisson_model <- tryCatch({
    glm(Count ~ Condition + Div, data = te_data, family = poisson)
  }, error = function(e) return(NULL))
  
  # Fit Negative Binomial Model
  nb_model <- tryCatch({
    glm.nb(Count ~ Condition + Div, data = te_data, control = glm.control(maxit = 25))
  }, error = function(e) return(NULL))
  
  # Fit Zero-Inflated Poisson Model
  zip_model <- tryCatch({
    zeroinfl(Count ~ Condition + Div | 1, data = te_data, dist = "poisson")
  }, error = function(e) return(NULL))
  
  # Fit Zero-Inflated Negative Binomial Model
  zinb_model <- tryCatch({
    zeroinfl(Count ~ Condition + Div | 1, data = te_data, dist = "negbin")
  }, error = function(e) return(NULL))
  
  # Store Results for Each Model
  aic_values <- c()
  p_values <- c()
  div_estimates <- c()
  overdispersion_values <- c()
  model_types <- c()
  shift_directions <- c()
  significance_stars <- c()
  
  # Check and store results for Poisson model
  if (!is.null(poisson_model)) {
    aic_values <- c(aic_values, AIC(poisson_model))
    model_types <- c(model_types, "Poisson")
    coef_summary <- summary(poisson_model)$coefficients
    if ("ConditionHeat_Stress" %in% rownames(coef_summary)) {
      p_values <- c(p_values, coef_summary["ConditionHeat_Stress", "Pr(>|z|)"])
      div_estimates <- c(div_estimates, coef_summary["Div", "Estimate"])
    } else {
      p_values <- c(p_values, NA)
      div_estimates <- c(div_estimates, NA)
    }
    overdispersion_values <- c(overdispersion_values, NA)  # Poisson has no overdispersion measure
    shift_directions <- c(shift_directions, ifelse(p_values[length(p_values)] < 0.05, "Shifted", "No Shift"))
    significance_stars <- c(significance_stars, ifelse(p_values[length(p_values)] < 0.001, "***",
                                                       ifelse(p_values[length(p_values)] < 0.01, "**",
                                                              ifelse(p_values[length(p_values)] < 0.05, "*", "ns"))))
  }
  
  # Check and store results for Negative Binomial model
  if (!is.null(nb_model)) {
    aic_values <- c(aic_values, AIC(nb_model))
    model_types <- c(model_types, "Negative Binomial")
    coef_summary <- summary(nb_model)$coefficients
    if ("ConditionHeat_Stress" %in% rownames(coef_summary)) {
      p_values <- c(p_values, coef_summary["ConditionHeat_Stress", "Pr(>|z|)"])
      div_estimates <- c(div_estimates, coef_summary["Div", "Estimate"])
    } else {
      p_values <- c(p_values, NA)
      div_estimates <- c(div_estimates, NA)
    }
    overdispersion_values <- c(overdispersion_values, overdispersion_test(nb_model))  # Calculate overdispersion
    shift_directions <- c(shift_directions, ifelse(p_values[length(p_values)] < 0.05, "Shifted", "No Shift"))
    significance_stars <- c(significance_stars, ifelse(p_values[length(p_values)] < 0.001, "***",
                                                       ifelse(p_values[length(p_values)] < 0.01, "**",
                                                              ifelse(p_values[length(p_values)] < 0.05, "*", "ns"))))
  }
  
  # Check and store results for Zero-Inflated Poisson model
  if (!is.null(zip_model)) {
    aic_values <- c(aic_values, AIC(zip_model))
    model_types <- c(model_types, "Zero-Inflated Poisson")
    coef_summary <- summary(zip_model)$coefficients
    if ("ConditionHeat_Stress" %in% rownames(coef_summary)) {
      p_values <- c(p_values, coef_summary["ConditionHeat_Stress", "Pr(>|z|)"])
      div_estimates <- c(div_estimates, coef_summary["Div", "Estimate"])
    } else {
      p_values <- c(p_values, NA)
      div_estimates <- c(div_estimates, NA)
    }
    overdispersion_values <- c(overdispersion_values, NA)  # Zero-Inflated Poisson has no overdispersion measure
    shift_directions <- c(shift_directions, ifelse(p_values[length(p_values)] < 0.05, "Shifted", "No Shift"))
    significance_stars <- c(significance_stars, ifelse(p_values[length(p_values)] < 0.001, "***",
                                                       ifelse(p_values[length(p_values)] < 0.01, "**",
                                                              ifelse(p_values[length(p_values)] < 0.05, "*", "ns"))))
  }
  
  # Check and store results for Zero-Inflated Negative Binomial model
  if (!is.null(zinb_model)) {
    aic_values <- c(aic_values, AIC(zinb_model))
    model_types <- c(model_types, "Zero-Inflated Negative Binomial")
    coef_summary <- summary(zinb_model)$coefficients
    if ("ConditionHeat_Stress" %in% rownames(coef_summary)) {
      p_values <- c(p_values, coef_summary["ConditionHeat_Stress", "Pr(>|z|)"])
      div_estimates <- c(div_estimates, coef_summary["Div", "Estimate"])
    } else {
      p_values <- c(p_values, NA)
      div_estimates <- c(div_estimates, NA)
    }
    overdispersion_values <- c(overdispersion_values, overdispersion_test(zinb_model))  # Calculate overdispersion
    shift_directions <- c(shift_directions, ifelse(p_values[length(p_values)] < 0.05, "Shifted", "No Shift"))
    significance_stars <- c(significance_stars, ifelse(p_values[length(p_values)] < 0.001, "***",
                                                       ifelse(p_values[length(p_values)] < 0.01, "**",
                                                              ifelse(p_values[length(p_values)] < 0.05, "*", "ns"))))
  }
  
  # Summarize Model Results in a DataFrame
  shift_results <- bind_rows(shift_results, tibble(
    TE_Family = te,
    Model_Type = model_types,
    AIC = aic_values,
    P_Value = p_values,
    Div_Estimate = div_estimates,
    Overdispersion = overdispersion_values,
    Shift_Direction = shift_directions,
    Significance = significance_stars
  ))
}

# Step 9: Save Model Summary
write.csv(shift_results, "TE_Family_Shift_Results_ovaries.csv", row.names = FALSE)

cat("\n✅ Analysis complete! Results saved to 'TE_Family_Shift_Results.csv'.\n")

# Plotting
# Load necessary libraries
library(tidyverse)
library(MASS)  # For Negative Binomial Model (glm.nb)
library(pscl)  # For Zero-Inflated Model (zeroinfl)

# 🔹 Step 1: Read Data
control_data <- read.delim("FC_all.divsum.csv", sep="", header=FALSE, check.names=FALSE, stringsAsFactors=FALSE)
temperature_data <- read.delim("FT_all.divsum.csv", sep="", header=FALSE, check.names=FALSE, stringsAsFactors=FALSE)

# 🔹 Step 2: Assign Column Names
colnames(control_data) <- as.character(control_data[1, ])
colnames(temperature_data) <- as.character(temperature_data[1, ])

# Remove the first row (now used as column names)
control_data <- control_data[-1, ]
temperature_data <- temperature_data[-1, ]

# Reset row numbers
rownames(control_data) <- NULL
rownames(temperature_data) <- NULL

# 🔹 Step 3: Convert Columns to Numeric & Remove Negatives
convert_numeric <- function(df) {
  df %>%
    mutate(across(where(is.character), ~ suppressWarnings(as.numeric(.)))) %>%
    replace(is.na(.), 0)  # Replace NA with 0
}

control_data <- convert_numeric(control_data)
temperature_data <- convert_numeric(temperature_data)

# Remove any negative values
control_data[control_data < 0] <- 0
temperature_data[temperature_data < 0] <- 0

cat("✅ No negative values remain.\n")

# 🔹 Step 4: Add Condition Column
control_data$Condition <- "Control"
temperature_data$Condition <- "Heat_Stress"

# 🔹 Step 5: Merge Data
merged_data <- bind_rows(control_data, temperature_data)

# 🔹 Step 6: Reshape Data to Long Format
data_long <- merged_data %>%
  pivot_longer(cols = -c(Div, Condition), names_to = "Repeat_Family", values_to = "Count")

# 🔹 Step 7: Extract TE Family Categories (DNA, LINE, LTR, SINE, RC)
data_long <- data_long %>%
  mutate(TE_Family = case_when(
    str_detect(Repeat_Family, "DNA")  ~ "DNA",
    str_detect(Repeat_Family, "LINE") ~ "LINE",
    str_detect(Repeat_Family, "LTR")  ~ "LTR",
    str_detect(Repeat_Family, "SINE") ~ "SINE",
    str_detect(Repeat_Family, "RC")   ~ "RC",
    str_detect(Repeat_Family, "Satellite")   ~ "Satellite",
    TRUE ~ "Other"
  ))

# 🔹 Step 8: Ensure Factors Are Correct
data_long$Condition <- as.factor(data_long$Condition)
data_long$TE_Family <- as.factor(data_long$TE_Family)

# 🔹 Step 9: Remove TE Families with Too Few Nonzero Counts
te_counts <- data_long %>%
  group_by(TE_Family) %>%
  summarise(NonZero_Counts = sum(Count > 0), Total_Counts = n())

# Keep only TE families with at least 200 nonzero counts
valid_te_families <- te_counts %>% filter(NonZero_Counts >= 200) %>% pull(TE_Family)
filtered_data <- data_long %>% filter(TE_Family %in% valid_te_families)

# 🔹 Step 10: Run Separate Models for Each TE Family
results <- list()
for (te in unique(filtered_data$TE_Family)) {
  
  cat("\n📌 Testing TE Family:", te, "...\n")
  
  te_data <- filtered_data %>% filter(TE_Family == te)
  
  # Fit Negative Binomial Model
  nb_model <- tryCatch({
    glm.nb(Count ~ Condition + Div, data = te_data, control = glm.control(maxit = 25))
  }, error = function(e) return(NULL))
  
  if (is.null(nb_model)) {
    cat("⚠️ Negative Binomial failed for", te, ". Trying Zero-Inflated Model...\n")
    
    zinb_model <- tryCatch({
      zeroinfl(Count ~ Condition + Div, data = te_data, dist = "negbin")
    }, error = function(e) return(NULL))
    
    if (!is.null(zinb_model)) {
      results[[te]] <- list(model = zinb_model, type = "Zero-Inflated NB")
    } else {
      cat("❌ Failed to fit any model for", te, "\n")
      next
    }
    
  } else {
    results[[te]] <- list(model = nb_model, type = "Negative Binomial")
  }
}

# 🔹 Step 11: Extract and Summarize P-Values for Condition Effect
summary_results <- tibble(
  TE_Family = names(results),
  Model_Type = sapply(results, function(x) x$type),
  P_Value = sapply(results, function(x) {
    coef_summary <- summary(x$model)$coefficients
    if ("ConditionHeat_Stress" %in% rownames(coef_summary)) {
      return(coef_summary["ConditionHeat_Stress", "Pr(>|z|)"])
    } else {
      return(NA)
    }
  })
)

# 🔹 Step 12: Sort and Display Results
summary_results <- summary_results %>% arrange(P_Value)

print(summary_results)

cat("\n✅ Analysis complete! Check p-values to see which TE families shift under heat stress.\n")


# check direction of shift in TE age
# Extract the direction of Kimura shift for each TE Family
direction_results <- tibble(
  TE_Family = names(results),
  Model_Type = sapply(results, function(x) x$type),
  Div_Estimate = sapply(results, function(x) {
    coef_summary <- summary(x$model)$coefficients
    if ("Div" %in% rownames(coef_summary)) {
      return(coef_summary["Div", "Estimate"])
    } else {
      return(NA)
    }
  }),
  P_Value_Div = sapply(results, function(x) {
    coef_summary <- summary(x$model)$coefficients
    if ("Div" %in% rownames(coef_summary)) {
      return(coef_summary["Div", "Pr(>|z|)"])
    } else {
      return(NA)
    }
  })
)

# Sort by significance (p-value) and print
direction_results <- direction_results %>% arrange(P_Value_Div)
print(direction_results)

# visual representation of shift in TE age
library(ggplot2)

ggplot(data_long, aes(x = Div, y = Count, color = Condition)) +
  geom_point(alpha=0.5) +
  geom_smooth(method = "loess") +
  facet_wrap(~ TE_Family, scales = "free_y") +
  theme_minimal() +
  scale_color_manual(values = c("Control" = "blue", "Heat_Stress" = "orange")) +  # Custom colors
  coord_cartesian(xlim = c(0, 30)) +  # Limit x-axis to 40
  labs(title = "TE Abundance vs. Kimura Divergence under Heat Stress",
       x = "Kimura Divergence (Div)",
       y = "TE Count",
       color = "Condition")

# lets add the negative binomial results
# Extract p-values and estimates for Div from models
div_significance <- tibble(
  TE_Family = names(results),
  Model_Type = sapply(results, function(x) x$type),
  Div_Estimate = sapply(results, function(x) {
    coef_summary <- summary(x$model)$coefficients
    if ("Div" %in% rownames(coef_summary)) {
      return(coef_summary["Div", "Estimate"])
    } else {
      return(NA)
    }
  }),
  P_Value_Div = sapply(results, function(x) {
    coef_summary <- summary(x$model)$coefficients
    if ("Div" %in% rownames(coef_summary)) {
      return(coef_summary["Div", "Pr(>|z|)"])
    } else {
      return(NA)
    }
  })
)

# Format p-values for display (convert to stars)
div_significance <- div_significance %>%
  mutate(Significance = case_when(
    P_Value_Div < 0.001 ~ "***",
    P_Value_Div < 0.01  ~ "**",
    P_Value_Div < 0.05  ~ "*",
    TRUE ~ "ns"
  ))

print(div_significance)  # Check results


# lets add this to our plot
library(ggplot2)

# Merge significance results into main dataset
data_long <- left_join(data_long, div_significance, by = "TE_Family")

# Create the plot with significance annotations
ggplot(data_long, aes(x = Div, y = Count, color = Condition)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  facet_wrap(~ TE_Family, scales = "free_y") +
  theme_minimal() +
  scale_color_manual(values = c("Control" = "blue", "Heat_Stress" = "red")) +  # Custom colors
  coord_cartesian(xlim = c(0, 40)) +  # Limit x-axis to 40
  labs(title = "TE Abundance vs. Kimura Divergence under Heat Stress",
       x = "Kimura Divergence (Div)",
       y = "TE Count",
       color = "Condition") +
  geom_text(data = div_significance, aes(x = 35, y = max(data_long$Count, na.rm = TRUE), 
                                         label = paste0("Div: ", Significance)), inherit.aes = FALSE)



# add direction of shift
# Create Younger/Older labels based on Div_Estimate
div_significance <- div_significance %>%
  mutate(Shift = case_when(
    P_Value_Div < 0.05 & Div_Estimate < 0 ~ "Younger",
    P_Value_Div < 0.05 & Div_Estimate > 0 ~ "Older",
    TRUE ~ "No Significant Shift"
  )) %>%
  mutate(Plot_Label = paste0(Shift, " (", Significance, ")"))  # Combine for plot

# update the plot
library(ggplot2)

# Merge significance results into main dataset for plotting
data_long <- left_join(data_long, div_significance, by = "TE_Family")

# Create the improved plot
ggplot(data_long, aes(x = Div, y = Count, color = Condition)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "loess") +
  facet_wrap(~ TE_Family, scales = "free_y") +
  theme_classic() +
  scale_color_manual(values = c("Control" = "blue", "Heat_Stress" = "orange")) +  # Custom colors
  coord_cartesian(xlim = c(0, 30), ylim = c(0, 100000000)) +  # Limit x-axis to 40
  labs(title = "Kimura substitution under Thermal Stress - Ovaries",
       x = "Kimura Substitution level",
       y = "TE abundance",
       color = "Condition") 
  # Add text annotation for Younger/Older shift at the top of each facet
  #geom_text(data = div_significance, 
   #         aes(x = 35, y = max(data_long$Count, na.rm = TRUE) * 0.9, label = Plot_Label), 
    #        inherit.aes = FALSE, color = "black", fontface = "bold", size = 5)

# don't plot other group TEs
ggplot(data_long %>% filter(TE_Family != "Other"), aes(x = Div, y = Count, color = Condition)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "loess") +
  facet_wrap(~ TE_Family, scales = "free_y") +
  theme_classic() +
  scale_color_manual(values = c("Control" = "blue", "Heat_Stress" = "orange")) +  # Custom colors
  coord_cartesian(xlim = c(0, 30), ylim = c(0, 100000000)) +  # Limit x-axis to 40
  labs(title = "Kimura substitution under Thermal Stress - Ovaries",
       x = "Kimura Substitution level",
       y = "TE abundance",
       color = "Condition") +
  theme(
    text = element_text(face = "bold", color = "black"),  # Make all text bold and black
    strip.text = element_text(face = "bold", color = "black"),  # Make facet labels bold and black
    axis.title = element_text(face = "bold", color = "black"),  # Make axis titles bold and black
    axis.text = element_text(face = "bold", color = "black")  # Make axis tick labels bold and black
  )


