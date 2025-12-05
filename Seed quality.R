# Seed Quality BLUPs

#Loading packages

library(lme4)
library(emmeans)
library(multcomp)
library(multcompView)

options(contrasts = c("contr.sum", "contr.sum"))

# Loading Data

seeds <- read.csv2("MP258_seed_quality.csv", stringsAsFactors = TRUE)


# Calculating the seed quality index 

# Defining good and bad traits:
good_traits <- c("Oil..DM", "Protein..DM", "C.18.1", "Sulphur..DM")
bad_traits  <- c("Glucosinolates.µmol.g.seed",
                 "Moisture..",
                 "C.18.3",
                 "C.22.1")

# Standardize the traits to make em comparable
all_traits <- c(good_traits, bad_traits)
z_all <- scale(seeds[, all_traits])
z_all <- as.data.frame(z_all)

# good traits positive, bad traits negative
sqi_per_row <- rowMeans(cbind(
  z_all[, good_traits],
  -z_all[, bad_traits]
))

seeds$SQI <- as.numeric(sqi_per_row)



# BLUP modell

# Genotype and Rep as random effects
m_blup <- lmer(SQI ~ 1 + (1 | Genotyp) + (1 | Rep), data = seeds)

# Extract the BLUPs
blup_values <- ranef(m_blup)$Genotyp
blup_df <- data.frame(
  Genotyp = rownames(blup_values),
  BLUP = blup_values[, 1]
)

# Calculating the grand mean
grand_mean <- fixef(m_blup)[1]

# EBLUP = BLUP + Grand Mean
blup_df$EBLUP <- blup_df$BLUP + grand_mean

# Ranking
blup_rank <- blup_df[order(-blup_df$EBLUP), ]
print(blup_rank)

# Highlighting our accessions
highlight <- c("ASSYST120", "ASSYST145", "ASSYST210",
               "ASSYST271", "ASSYST419", "ASSYST424")

colors_blup <- ifelse(blup_rank$Genotyp %in% highlight,
                      "darkgreen", "lightblue")

# BLUP-Plot
par(mar = c(10, 5, 4, 2))
barplot(
  blup_rank$EBLUP,
  names.arg = blup_rank$Genotyp,
  las = 2,
  col = colors_blup,
  cex.names = 0.8,
  ylab = "Seed Quality (BLUP)",
  main = "BLUP-based Seed Quality Ranking"
)



# Emmeans just in case


seeds$G <- as.factor(seeds$Genotyp) 
seeds$R <- as.factor(seeds$Rep)

m <- lmer(SQI ~ G + (1|R), data = seeds)

anova(m)

# Estimated marginal means per genotype
emm <- emmeans(m, ~ G)

# Tukey Test
pairs(emm, adjust = "tukey")

# Compact Letter Display
cld_tab <- cld(emm, Letters = letters, alpha = 0.1)
cld_tab <- as.data.frame(cld_tab)

# whoopsie
colnames(cld_tab)[colnames(cld_tab) == "G"] <- "Genotyp"

# sort after quality
cld_tab <- cld_tab[order(-cld_tab$emmean), ]

# change colors
cld_tab$col <- ifelse(cld_tab$Genotyp %in% highlight,
                      "darkgreen", "lightblue")

cld_tab

# Plotting this...

par(mar = c(10, 5, 4, 2))

bp <- barplot(
  cld_tab$emmean,
  names.arg = cld_tab$Genotyp,
  las = 2,
  col = cld_tab$col,
  ylab = "Estimated SQI",
  main = "Seed Quality Index with Significance Groups",
  cex.names = 0.8,
  ylim = c(min(cld_tab$emmean) - 0.1,
           max(cld_tab$emmean) + 0.15)
)

# I also want the groups
text(
  x = bp,
  y = cld_tab$emmean,
  labels = gsub(" ", "", cld_tab$.group),
  pos = ifelse(cld_tab$emmean < 0, 1, 3),
  cex = 1.2,
  font = 2
)

mtext("Genotype", side = 1, line = 8)
