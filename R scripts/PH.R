
# Packages
library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(multcomp)
library(multcompView)
library(emmeans)

# Loading Data
df <- read.csv2("PH.csv", stringsAsFactors = TRUE)

#Wide to Long
df_long <- df %>%
  pivot_longer(
    cols = c(PH1, PH2),  # column names
    names_to = "Rep",               # new column for rep
    values_to = "PH"               # new column for trait value
  )

df_long$Rep <- as.factor(df_long$Rep)

head(df_long)

# Looks good

model_fixed <- lmer(PH ~ Genotype + (1|Rep), data = df_long)

anova(model_fixed)

emm <- emmeans(model_fixed, ~ Genotype)


pairs(emm, adjust = "tukey")


cld(emm,
    alpha = 0.05,
    Letters = letters,
    adjust = "tukey"  
)

# genotype is significant

# BLUP and EBLUP calculations

model_blup <- lmer(PH ~ 1 + (1|Genotype) + (1|Rep), data = df_long)

bl <- ranef(model_blup)$Genotype
bl <- data.frame(
  Genotype = rownames(bl),
  BLUP = bl[,1]
)

grand_mean <- fixef(model_blup)[1]
bl$EBLUP <- bl$BLUP + grand_mean

print(bl)