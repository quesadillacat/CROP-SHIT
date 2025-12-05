# Flowering Time

# Ok... So... this one is... hierarchical???????? Because genotype is like
# nested in type? environment is random because we are looking at the genotypes
# and not at the environments.. right?! I mean, we want to know which genotype
# is better and not which location?


library(tidyr)
library(dplyr)
library(lme4)
library(emmeans)
library(multcomp)
library(multcompView)

# Load se data

ft <- read.csv2("MK130_OPTIONAL_flowering_time_data.csv", stringsAsFactors = FALSE)

# Genotype factor defined
ft$Genotype <- as.factor(ft$Accession_Number)

# Type-factor defined
ft$Type <- as.factor(ft$type)

# Convert to proper long format

ft_long <- ft %>%
  pivot_longer(
    cols = matches("GI|GG|RH|CQ|TE"),
    names_to = "Environment",
    values_to = "FT"
  ) %>%
  filter(!is.na(FT))

# Every environment column has location and year
# ...
# She could have at least... separated them before asking us to analyze that shit


ft_long$Environment <- as.factor(ft_long$Environment)
ft_long$Location    <- as.factor(substr(ft_long$Environment, 1, 2))
ft_long$Year        <- as.factor(substr(ft_long$Environment, 3, 6))

# Mixed Model: Type + Genotype NESTED in (Type) + random Environment

m2 <- lmer(FT ~ Type + Genotype:Type + (1 | Location:Year),
           data = ft_long)

anova(m2)


# Estimated Means per TYPE

emm_type <- emmeans(m2, ~ Type)
emm_type
cld_type <- cld(emm_type, Letters = letters, alpha = 0.1)
cld_type

#  Estimated Means per Genotype WITHIN each Type

emm_geno <- emmeans(m2, ~ Genotype | Type)
cld_geno <- cld(emm_geno, Letters = letters, alpha = 0.1)

# cld displayy

cld_geno <- as.data.frame(cld_geno)
cld_geno$.group <- gsub(" ", "", cld_geno$.group)

# Sort within each type by flowering time

cld_geno <- cld_geno %>%
  group_by(Type) %>%
  arrange(emmean, .by_group = TRUE)

cld_geno

# Ok so ig Genotype doesn't matter for flowering time... not relevant
# for our choice?? Not relevant for PCA or multi-trait-score

# I am VERY confused about this one tho. did i do it correctly? i mean i did
# do the means for the genotypes within each of the types. so it should be good??

# Means for PCA???
FT_means <- ft_long %>%
  group_by(Genotype) %>%
  summarise(FT_mean = mean(FT, na.rm = TRUE))

FT_means
