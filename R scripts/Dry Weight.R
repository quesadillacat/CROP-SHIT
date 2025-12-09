
# Packages
library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(multcompView)

# Loading Data
df <- read.csv2("DW.csv", stringsAsFactors = TRUE)

#Wide to Long
df_long <- df %>%
  pivot_longer(
    cols = c(DW1, DW2),  # column names
    names_to = "Rep",               # new column for rep
    values_to = "DW"               # new column for trait value
  )

df_long$Rep <- as.factor(df_long$Rep)

head(df_long)

# Looks good

# Does Genotype influence CCM?

model_dw_fixed <- lm(DW ~ Genotype, data = df_long)

anova(model_dw_fixed)

emm_dw <- emmeans(model_dw_fixed, ~ Genotype)

pairs(emm_dw, adjust = "tukey")

# cld

cld_dw <- multcomp::cld(
  emm_dw,
  adjust = "tukey",
  Letters = letters
)

print(cld_dw)

# It does, big time!!!!

# BLUP and EBLUP calculations

model_blup <- lmer(DW ~ 1 + (1|Genotype) + (1|Rep), data = df_long)

bl <- ranef(model_blup)$Genotype
bl <- data.frame(
  Genotype = rownames(bl),
  BLUP = bl[,1]
)

grand_mean <- fixef(model_blup)[1]
bl$EBLUP <- bl$BLUP + grand_mean

print(bl)