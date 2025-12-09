
# Packages
library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)

# Loading Data
df <- read.csv2("HGW.csv", stringsAsFactors = TRUE)

#Wide to Long
df_long <- df %>%
  pivot_longer(
    cols = c(HGW_1, HGW_2, HGW_3, HGW_4, HGW_5, HGW_6),  # column names
    names_to = "Rep",               # new column for rep
    values_to = "HGW"               # new column for trait value
  )

df_long$Rep <- as.factor(df_long$Rep)

head(df_long)

# Looks good

# Does Genotype influence CCM?

model_fixed <- lmer(HGW ~ Genotype + (1|Rep), data = df_long)

anova(model_fixed)

# yaaay it does

# BLUP and EBLUP calculations

model_blup <- lmer(HGW ~ 1 + (1|Genotype) + (1|Rep), data = df_long)

bl <- ranef(model_blup)$Genotype
bl <- data.frame(
  Genotype = rownames(bl),
  BLUP = bl[,1]
)

grand_mean <- fixef(model_blup)[1]
bl$EBLUP <- bl$BLUP + grand_mean

print(bl)


