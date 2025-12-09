# Load packages
library(tidyr) # Load data into long format
library(dplyr) 
library(lme4)  #Linear mixed model
library(lmerTest)  #p-values for lmms
library(multcomp)
library(multcompView)
library(emmeans)

# Load Data
df <- read.csv2("BBCH.csv", stringsAsFactors = TRUE) #cvs.2 uses ; as separator

# Convert table to long format
df_long <- df %>%
  pivot_longer(                  #pivot_longer makes one column out of BBCH1 and BBCH2
    cols = c(BBCH1, BBCH2),   
    names_to = "Rep",
    values_to = "BBCH"
  )

# Define Factors
df_long$Rep      <- as.factor(df_long$Rep)
df_long$Genotype <- as.factor(df_long$Genotype)
df_long$Type     <- as.factor(df_long$Type)


# Are genotypes significant WITHOUT type effect?

model_geno_fixed <- lm(
  BBCH ~ Type + Genotype,
  data = df_long
)

anova(model_geno_fixed)

# Tukey
emmeans(model_geno_fixed, pairwise ~ Genotype, adjust = "tukey")

# cld
emm_global <- emmeans(model_geno_fixed, ~ Genotype)

cld_global <- multcomp::cld(
  emm_global,
  adjust = "tukey",
  Letters = letters
)

cld_global

#Yup, they are.That's why this belongs to the multi-trait-score!

# Is Genotype significant within type? 
df_long$Genotype_in_Type <- interaction(df_long$Type, df_long$Genotype)

model_nested <- lm(BBCH ~ Type + Type:Genotype, data = df_long)


anova(model_nested)
# Tukey
emmeans(model_nested, pairwise ~ Genotype | Type, adjust = "tukey")

# cld
emm_nested <- emmeans(model_nested, ~ Genotype | Type)

cld_nested <- multcomp::cld(
  emm_nested,
  adjust = "tukey",
  Letters = letters
)

cld_nested

# Well, barely.

# Calculate BLUP:

model_blup_type <- lmer(
  BBCH ~ Type + (1|Genotype) + (1|Rep),
  data = df_long
)

# EXTRACT GENOTYPE BLUPS
bl <- ranef(model_blup_type)$Genotype
bl <- data.frame(
  Genotype = rownames(bl),
  BLUP = bl[,1]
)

# GRAND MEAN 
grand_mean <- fixef(model_blup_type)[["(Intercept)"]]

# EBLUP
bl$EBLUP <- bl$BLUP + grand_mean

# OUTPUT 
print(bl)


