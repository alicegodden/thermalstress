# Title: GLM testing on Retroseq data
# Author: Dr. Alice M. Godden

# Author : Sara Irish
library(dplyr)

library(tidyr)

library(lubridate)

library(survival)

library(survminer)

library(ggplot2)

library(dabestr)

library(lme4)

library(coxme)

library(lmerTest)

library(glmmTMB)

library(emmeans)

library(DHARMa)

library(car)

# Human Methylcellulose data

# read in your data
df <- read.csv("zebrafish_TE_counts_heat_stress_te.csv")
df
head (df)
str(df)


#repro_long is the name of the new longform dataset

repro_long <- df %>% 
  
  pivot_longer(
    
    cols = c(`Total`,`DNA`, `LINE`, `SINE`, `LTR`, `SATELLITE`, `RC`), #list all TE columns in dataset here RC for fish
    
    names_to = "TE_type", #this creates a column named TE_type, containing Total and TE by family
    
    values_to = "value") #this creates a column named value which contains the offspring numbers produced each day



head(repro_long)




# plot the TE data
repro_plot <- ggplot(data=repro_long, aes(x=TE_type, y=value, group=Treatment, color=Treatment)) +
  geom_jitter(alpha = 0.2, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", size = 1.2, width=0.0, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="point", size = 3, position = position_dodge(0.5)) +
  stat_summary(fun.data="mean_cl_boot", geom="line",  size=1.2, position = position_dodge(0.5)) +
  theme_classic() +
  labs(y="TE count", x="TE", title ="Zebrafish- TE counts- unique to treatment group") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.position = c(0.8, 0.9)
  ) +
  scale_color_manual(values = c(
    "Testes_Ctrl" = "blue",
    "Testes_Temp" = "red"
  ))

# generate plot
repro_plot

# Save the plot to a file
ggsave(
  filename = "zfish_TE_count_plot_HS_te.png",   # The name of your file and its extension
  plot = repro_plot,            # The ggplot object you want to save
  width = 12,                    # Width in inches (adjust as needed)
  height = 5,                   # Height in inches (adjust as needed)
  dpi = 600                     # Resolution (dots per inch) for raster images
)

# now make a histogram
hist <- hist(repro_long$value)
ggsave(
  filename = "zfish_TE_count_hist_HS_ov.png",   # The name of your file and its extension
  plot = hist,            # The ggplot object you want to save
  width = 8,                    # Width in inches (adjust as needed)
  height = 6,                   # Height in inches (adjust as needed)
  dpi = 600                     # Resolution (dots per inch) for raster images
)


#models
repro_long$Time <- as.factor(repro_long$Treatment)
##Then, change the variables names to numbers
levels(repro_long$Treatment) <- list('1' = 'Ovaries_Ctrl', '2' = 'Ovaries_Temp') # for MC-'1' = 'raw', '2' = 'central', '3' ='outer'
levels(repro_long$Treatment)

# running poisson model
rep_m1 <- glmer(value ~ Treatment + (1 | MaleID), family = poisson, data = repro_long)
summary(rep_m1)


# try a negative binomial if there is overdispersion

# Fit a Negative Binomial GLMM
rep_m_negbin <- glmer.nb(value ~ Treatment + (1 | MaleID), data = repro_long)
summary(rep_m_negbin)

# analyse if there is a significant change in TE by family/group
df$Treatment <- as.factor(df$Treatment)
levels(df$Treatment) <- list('1' = 'Testes_Ctrl', '2' = 'Testes_Temp') #, '3' ='outer') 1 is control grouup/ref group
levels(df$Treatment)


# running poisson model
rep_m1 <- glmer(Total ~ Treatment + (1 | MaleID), family = poisson, data = df)
summary(rep_m1)

# Get the residual deviance and degrees of freedom from the model summary
library(performance)
check_overdispersion(rep_m1)
check_overdispersion(rep_m_negbin_SATELLITTE)

# try a negative binomial if there is overdispersion
rep_m_negbin <- glmer.nb(Total ~ Treatment + (1 | MaleID), data = df)
summary(rep_m_negbin)

# Fit a Negative Binomial GLMM
rep_m_negbin_total <- glmer.nb(Total ~ Treatment + (1 | MaleID), data = df)
summary(rep_m_negbin_total)

rep_m_negbin_DNA <- glmer.nb(DNA ~ Treatment + (1 | MaleID), data = df)
summary(rep_m_negbin_DNA)

rep_m_negbin_LINE <- glmer.nb(LINE ~ Treatment + (1 | MaleID), data = df)
summary(rep_m_negbin_LINE)

rep_m_negbin_SINE <- glmer.nb(SINE ~ Treatment + (1 | MaleID), data = df)
summary(rep_m_negbin_SINE)

rep_m_negbin_LTR <- glmer.nb(LTR ~ Treatment + (1 | MaleID), data = df)
summary(rep_m_negbin_LTR)

rep_m_negbin_SATELLITTE <- glmer.nb(SATELLITE ~ Treatment + (1 | MaleID), data = df)
summary(rep_m_negbin_SATELLITTE)


rep_m_negbin_RC <- glmer.nb(RC ~ Treatment + (1 | MaleID), data = df)
summary(rep_m_negbin_RC)
# fit poisson
# Change these Negative Binomial GLMMs to Poisson GLMMs

# Fit a Poisson GLMM for Total
check_overdispersion(rep_m_poisson_SATELLITTE)
rep_m_poisson_total <- glmer(Total ~ Treatment + (1 | MaleID), family = poisson, data = df)
summary(rep_m_poisson_total)

# Fit a Poisson GLMM for DNA
rep_m_poisson_DNA <- glmer(DNA ~ Treatment + (1 | MaleID), family = poisson, data = df)
summary(rep_m_poisson_DNA)

# Fit a Poisson GLMM for LINE
rep_m_poisson_LINE <- glmer(LINE ~ Treatment + (1 | MaleID), family = poisson, data = df)
summary(rep_m_poisson_LINE)

# Fit a Poisson GLMM for SINE
rep_m_poisson_SINE <- glmer(SINE ~ Treatment + (1 | MaleID), family = poisson, data = df)
summary(rep_m_poisson_SINE)

# Fit a Poisson GLMM for LTR
rep_m_poisson_LTR <- glmer(LTR ~ Treatment + (1 | MaleID), family = poisson, data = df)
summary(rep_m_poisson_LTR)

# Fit a Poisson GLMM for SATELLITTE
rep_m_poisson_SATELLITTE <- glmer(SATELLITE ~ Treatment + (1 | MaleID), family = poisson, data = df)
summary(rep_m_poisson_SATELLITTE)


rep_m_poisson_RC <- glmer(RC ~ Treatment + (1 | MaleID), family = poisson, data = df)
summary(rep_m_poisson_RC)

# running to see central vs outer 2
# --- 1. Set 'Time2' as the reference level for the 'Time' factor ---
# Assuming your Time levels are '1', '2', '3'
df$Time <- relevel(df$Time, ref = "1")

#then re-run the above

library(lme4)
summary(lmer(Total ~ Time + (1 | MaleID), data = df))
