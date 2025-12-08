
# Packages
library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)

# Loading Data
df <- read.csv2("CCM.csv", stringsAsFactors = TRUE)

#Wide to Long
df_long <- df %>%
  pivot_longer(
    cols = c(CCM_1, CCM_2, CCM_3, CCM_4, CCM_5, CCM_6),  # column names
    names_to = "Rep",               # new column for rep
    values_to = "CCM"               # new column for trait value
  )

df_long$Rep <- as.factor(df_long$Rep)

head(df_long)

# Looks good

# Does Genotype influence CCM?

# Use a fixed model
model_ccm_fixed <- lm(CCM ~ Genotype, data = df_long)

# ANOVA
anova(model_ccm_fixed)

# Tukey and emmeans
emm_ccm <- emmeans(model_ccm_fixed, ~ Genotype)

pairs(emm_ccm, adjust = "tukey")

# CLD
cld_ccm <- multcomp::cld(
  emm_ccm,
  adjust = "tukey",
  Letters = letters
)

print(cld_ccm)

# yaaay it does

# BLUP and EBLUP calculations

model_blup <- lmer(CCM ~ 1 + (1|Genotype) + (1|Rep), data = df_long)

bl <- ranef(model_blup)$Genotype
bl <- data.frame(
  Genotype = rownames(bl),
  BLUP = bl[,1]
)

grand_mean <- fixef(model_blup)[1]
bl$EBLUP <- bl$BLUP + grand_mean

print(bl)


