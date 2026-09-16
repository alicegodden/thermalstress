# Testing TE counts class unique counts to thermal groups

# =====================================================================
# Zebrafish de novo TE insertions (treatment-unique) under heat stress
# Revised Fig 3E/F + GLMM stats (Poisson / negative binomial per family)
# Handles testes (MaleID) or ovaries (FemaleID) input.
# Stats retained from S. Irish's original script.
# =====================================================================
library(dplyr); library(tidyr); library(ggplot2)
library(lme4); library(performance)

# ---- read BOTH organs and combine -----------------------------------
te <- read.csv("zebrafish_TE_counts_heat_stress_te.csv", check.names = TRUE)
ov <- read.csv("zebrafish_TE_counts_heat_stress_ov.csv", check.names = TRUE)

# the id column is MaleID in testes, FemaleID in ovaries; standardise to SampleID
id_col <- function(d) names(d)[grepl("ID$", names(d))][1]
names(te)[names(te) == id_col(te)] <- "SampleID"
names(ov)[names(ov) == id_col(ov)] <- "SampleID"
te$Organ <- "Testes"; ov$Organ <- "Ovaries"
df <- bind_rows(te, ov)
df$Cond <- sub(".*_", "", df$Treatment)        # Ctrl / Temp
fams <- c("Total","DNA","LINE","SINE","LTR","SATELLITE","RC")

# =====================================================================
# FIGURE (revised): TE count by family, testes vs ovaries
# =====================================================================
long <- df %>% pivot_longer(cols = all_of(fams),
                            names_to = "TE_type", values_to = "value") %>%
  mutate(TE_type = factor(TE_type, levels = fams),
         Group = paste(Organ, Cond))

# --- Option A: all four groups (keeps control baseline; recommended) ---
pA <- ggplot(long, aes(TE_type, value, group = Group, color = Group)) +
  geom_jitter(alpha = .2, position = position_jitterdodge(jitter.width = .2, dodge.width = .6)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", linewidth = 1, width = 0, position = position_dodge(.6)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 3, position = position_dodge(.6)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "line", linewidth = 1, position = position_dodge(.6)) +
  scale_color_manual(values = c("Testes Ctrl" = "#2c5f9e", "Testes Temp" = "#c85a28",
                                "Ovaries Ctrl" = "#5aa06a", "Ovaries Temp" = "#d98f52")) +
  labs(y = "TE count (treatment-unique insertions)", x = "TE family",
       title = "De novo TE insertions unique to each group") +
  theme_classic() +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text  = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11, face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave("fig3EF_TE_fourgroup.png", pA, width = 12, height = 5, dpi = 600)

# --- Option B: temperature-unique only, testes vs ovaries (reviewer's ask) ---
pB <- ggplot(filter(long, Cond == "Temp"),
             aes(TE_type, value, group = Organ, color = Organ)) +
  geom_jitter(alpha = .2, position = position_jitterdodge(jitter.width = .2, dodge.width = .5)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", linewidth = 1.1, width = 0, position = position_dodge(.5)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 3.2, position = position_dodge(.5)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "line", linewidth = 1.1, position = position_dodge(.5)) +
  scale_color_manual(values = c("Testes" = "#c85a28", "Ovaries" = "#d98f52")) +
  labs(y = "Temperature-unique TE count", x = "TE family",
       title = "Heat-mobilized (temperature-unique) TE insertions: testes vs ovaries") +
  theme_classic() +
  theme(plot.title = element_text(size = 15, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text  = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(face = "bold"))
ggsave("fig3EF_TE_temp_only.png", pB, width = 11, height = 5, dpi = 600)

# =====================================================================
# STATS — retained from original: per-family GLMM, Ctrl vs Temp,
# random intercept for sample (MaleID / FemaleID -> SampleID).
# Run separately per organ so the contrast is within-sex.
# =====================================================================
fit_family <- function(d, organ) {
  d <- d[d$Organ == organ, ]
  d$Treatment <- factor(d$Cond, levels = c("Ctrl","Temp"))  # Ctrl = reference
  cat("\n==================== ", organ, " ====================\n")
  for (fam in fams) {
    cat("\n----", fam, "----\n")
    f_pois <- as.formula(paste0(fam, " ~ Treatment + (1 | SampleID)"))
    m_p <- try(glmer(f_pois, family = poisson, data = d), silent = TRUE)
    if (!inherits(m_p, "try-error")) {
      od <- try(check_overdispersion(m_p), silent = TRUE)
      print(summary(m_p)$coefficients)
      if (!inherits(od,"try-error")) print(od)
      # negative binomial if overdispersed
      if (!inherits(od,"try-error") && !is.null(od$p_value) && od$p_value < 0.05) {
        cat("  -> overdispersed; fitting negative binomial\n")
        m_nb <- try(glmer.nb(f_pois, data = d), silent = TRUE)
        if (!inherits(m_nb,"try-error")) print(summary(m_nb)$coefficients)
      }
    } else cat("  model failed to converge\n")
  }
}
fit_family(df, "Testes")
fit_family(df, "Ovaries")

message("Done: fig3EF_TE_fourgroup.png, fig3EF_TE_temp_only.png + GLMM output above")
