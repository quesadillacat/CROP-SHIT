library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(multcompView)

df <- read.csv2("FM.csv", stringsAsFactors = TRUE)

df_long <- df %>%
  pivot_longer(
    cols = c(FM1, FM2),
    names_to = "Rep",
    values_to = "FM"
  )

df_long$Rep <- as.factor(df_long$Rep)

model_fixed <- lmer(FM ~ Genotype + (1|Rep), data = df_long)

anova(model_fixed)

emm <- emmeans(model_fixed, ~ Genotype)
pairs(emm, adjust = "tukey")

cld(emm,
    alpha = 0.05,
    Letters = letters,
    adjust = "tukey"
)

model_blup <- lmer(FM ~ 1 + (1|Genotype) + (1|Rep), data = df_long)

bl <- ranef(model_blup)$Genotype
bl <- data.frame(
  Genotype = rownames(bl),
  BLUP = bl[,1]
)

grand_mean <- fixef(model_blup)[1]
bl$EBLUP <- bl$BLUP + grand_mean

print(bl)
