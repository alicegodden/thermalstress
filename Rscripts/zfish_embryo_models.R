############################################################
# ZEBRAFISH FERTILISATION, BREEDING AND EMBRYO MORTALITY

# Author: Sara Irish, Alice Godden & Chelsea Drake
############################################################

############################################################
# 1. LOAD PACKAGES
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(ggbeeswarm)

############################################################
# 2. IMPORT DATA
############################################################

df <- read.csv("zebratemp.csv")

# Clean spelling of FemaleTemp
df <- df %>% rename(FemaleTemp = FeamleTemp)

# Create key proportions
df <- df %>%
  mutate(
    PropFert = Fertilised / TotalEggs,
    MortProp = Dead / Fertilised
  )

############################################################
# 3. CREATE CONSISTENT TREATMENT FACTOR
############################################################

df$Treatment <- factor(
  paste0(df$FemaleTemp, "F x ", df$MaleTemp, "M"),
  levels = c("28F x 28M", "28F x 34M", "34F x 28M", "34F x 34M")
)

############################################################
# 4. HELPER FUNCTION: Standardised barplot with SE bars
############################################################

plot_bar <- function(data, y, ylab, fill_by = "Treatment") {
  ggplot(data, aes_string(x = "Treatment", y = y, fill = fill_by)) +
    stat_summary(fun = mean, geom = "bar", colour = "black") +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    scale_fill_manual(values = c(
      "28F x 28M" = "cadetblue2",
      "28F x 34M" = "chocolate1",
      "34F x 28M" = "pink",
      "34F x 34M" = "darkseagreen2"
    )) +
    ylab(ylab) +
    xlab("Parental Temperature Treatment") +
    theme_classic() +
    theme(legend.position = "none")
}

############################################################
# 5. MODEL 1: FERTILISATION GLM
############################################################

fert_model <- glm(
  cbind(Fertilised, TotalEggs - Fertilised) ~ FemaleTemp * MaleTemp + Block,
  family = binomial,
  data = df
)
summary(fert_model)

# Diagnostic check
simulationOutput <- simulateResiduals(fert_model)
plot(simulationOutput)

# Plot
plot_bar(df, "PropFert", "Proportion fertilised (2 h)")

############################################################
# 6. MODEL 2: MORTALITY GLM
############################################################

mort_model <- glm(
  cbind(Dead, Fertilised - Dead) ~ FemaleTemp * MaleTemp + Block,
  family = binomial,
  data = df
)
summary(mort_model)

# Diagnostics
simulationOutput <- simulateResiduals(mort_model)
plot(simulationOutput)

# Plot
plot_bar(df, "MortProp", "Mortality Rate (24 h)")

############################################################
# 7. MODEL 3: TOTAL EGG PRODUCTION (QUASIPOISSON)
############################################################

egg_model <- glm(
  TotalEggs ~ FemaleTemp * MaleTemp + Block,
  family = quasipoisson,
  data = df
)
summary(egg_model)

# Diagnostics
simulationOutput <- simulateResiduals(egg_model)
plot(simulationOutput)

# Plot
plot_bar(df, "TotalEggs", "Mean Total Eggs Produced")

############################################################
# 8. SECOND DATASET: BREEDING SUCCESS
############################################################

df1 <- read.csv("zebrafishtotal.csv")
df1 <- df1 %>% rename(FemaleTemp = FeamleTemp)

df1$Bred <- as.numeric(df1$TotalEggs > 0)
df1$Treatment <- factor(
  paste0(df1$FemaleTemp, "F x ", df1$MaleTemp, "M"),
  levels = c("28F x 28M", "28F x 34M", "34F x 28M", "34F x 34M")
)

############################################################
# 9. MODEL 4: BREEDING PROBABILITY
############################################################

breed_model <- glm(
  Bred ~ FemaleTemp * MaleTemp + Block,
  family = binomial,
  data = df1
)
summary(breed_model)

# Diagnostics
simulationOutput <- simulateResiduals(breed_model)
plot(simulationOutput)

# Plot breeding success
plot_bar(df1, "Bred", "Proportion of pairs that bred")

############################################################
# 10. Optional beeswarm: raw mortality counts
############################################################

ggplot(df, aes(x = Treatment, y = Dead, colour = Treatment)) +
  geom_beeswarm(size = 2) +
  scale_colour_manual(values = c(
    "28F x 28M" = "cadetblue2",
    "28F x 34M" = "chocolate1",
    "34F x 28M" = "pink",
    "34F x 34M" = "darkseagreen2"
  )) +
  ylab("Number of Dead Embryos (24 h)") +
  xlab("Parental Temperature Treatment") +
  theme_classic() +
  theme(legend.position = "none")
