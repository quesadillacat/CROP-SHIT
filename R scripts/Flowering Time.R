library(lme4)
library(emmeans)
library(multcomp)
library(multcompView)
library(dplyr)
library(lmerTest)

ft <- read.csv2("MK130_OPTIONAL_flowering_time_data.csv", stringsAsFactors = FALSE)


# Define factors 
ft$Genotype <- as.factor(ft$Accession_Number) 
ft$Type <- as.factor(ft$type)

ft_long <- ft %>%
  pivot_longer(
    cols = matches("GI|GG|RH|CQ|TE"),
    names_to = "Environment",
    values_to = "FT"
  ) %>%
  filter(!is.na(FT))

# Factors after pivot
ft_long$Environment <- as.factor(ft_long$Environment)
ft_long$Location <- as.factor(substr(ft_long$Environment, 1, 2))
ft_long$Year <- as.factor(substr(ft_long$Environment, 3, 6))
# I'm dumb as fuck idk how to do loops......... sorry

# For EX -------------------------------------------------------------------
dat_EX <- subset(ft_long, Type == "EX")

# Modell: Genotype fix + Location:Year als random
m_EX <- lmer(FT ~ Genotype + (1 | Location:Year), data = dat_EX)

# ANOVA
anova(m_EX)

# EMMs
emm_EX <- emmeans(m_EX, ~ Genotype)

# Tukey
pairs(emm_EX, adjust = "tukey")

# cld
cld_EX <- cld(
  emm_EX,
  Letters = letters,
  adjust  = "tukey",
  alpha   = 0.05
)

cld_EX

# For SWE -----------------------------------------------------------------

dat_SWE <- subset(ft_long, Type == "SWE")

# Modell: Genotype fix + Location:Year als random
m_SWE <- lmer(FT ~ Genotype + (1 | Location:Year), data = dat_SWE)

# ANOVA
anova(m_SWE)

# EMMs
emm_SWE <- emmeans(m_SWE, ~ Genotype)

# Tukeyy
pairs(emm_SWE, adjust = "tukey")

# CLD
cld_SWE <- cld(
  emm_SWE,
  Letters = letters,
  adjust  = "tukey",
  alpha   = 0.05
)

cld_SWE

# For SOSR -----------------------------------------------------------------

dat_SOSR <- subset(ft_long, Type == "SOSR")

# Modell: Genotype fix + Location:Year als random
m_SOSR <- lmer(FT ~ Genotype + (1 | Location:Year), data = dat_SOSR)

# ANOVA
anova(m_SOSR)

# EMMs
emm_SOSR <- emmeans(m_SOSR, ~ Genotype)

# Tukey
pairs(emm_SOSR, adjust = "tukey")

# CLD
cld_SOSR <- cld(
  emm_SOSR,
  Letters = letters,
  adjust  = "tukey",
  alpha   = 0.05
)

cld_SOSR

# For WOSR ----------------------------------------------------------------

dat_WOSR <- subset(ft_long, Type == "WOSR")

# Modell: Genotype fix + Location:Year als random
m_WOSR <- lmer(FT ~ Genotype + (1 | Location:Year), data = dat_WOSR)

# ANOVA
anova(m_WOSR)

# EMMs
emm_WOSR <- emmeans(m_WOSR, ~ Genotype)

#Tukey
pairs(emm_WOSR, adjust = "tukey")

#Cld
cld_WOSR <- cld(
  emm_WOSR,
  Letters = letters,
  adjust  = "tukey",
  alpha   = 0.05
)

cld_WOSR
