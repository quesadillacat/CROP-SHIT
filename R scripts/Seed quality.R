# Seed Quality BLUPs – corrected version
library(lme4)
library(emmeans)
library(multcomp)
library(multcompView)

options(contrasts = c("contr.sum", "contr.sum"))

# Loading Data
seeds <- read.csv2("MP258_seed_quality.csv", stringsAsFactors = TRUE)

# Defining traits according to Prof Sarah

good_traits <- c(
  "Oil..DM",
  "Protein..DM",
  "C.18.1",
  "C.18.3"       
)

bad_traits <- c(
  "Glucosinolates.µmol.g.seed",
  "Moisture..",
  "C.22.1",
  "Sulphur..DM"   # 
)

all_traits <- c(good_traits, bad_traits)

# Check if all columns exist:
stopifnot(all(all_traits %in% names(seeds)))

# standardizing the traits

z_all <- scale(seeds[, all_traits])
z_all <- as.data.frame(z_all)

# Calculating the SQI

sqi_per_row <- rowMeans(cbind(
  z_all[, good_traits],     # good = positive
  -z_all[, bad_traits]      # bad = inverted
))

seeds$SQI <- as.numeric(sqi_per_row)

# Calculating and extracting BLUPs

m_blup <- lmer(SQI ~ 1 + (1 | Genotyp) + (1 | Rep),
               data = seeds)

blup_values <- ranef(m_blup)$Genotyp

blup_df <- data.frame(
  Genotyp = rownames(blup_values),
  BLUP = blup_values[,1]
)
# EBLUPs
grand_mean <- fixef(m_blup)[1]
blup_df$EBLUP <- blup_df$BLUP + grand_mean

# Ranking
blup_rank <- blup_df[order(-blup_df$EBLUP), ]
print(blup_rank)

# Plotting just in case.

highlight <- c("ASSYST120","ASSYST145","ASSYST210",
               "ASSYST271","ASSYST419","ASSYST424")

colors_blup <- ifelse(blup_rank$Genotyp %in% highlight,
                      "darkgreen", "lightblue")

par(mar = c(10,5,4,2))
barplot(
  blup_rank$EBLUP,
  names.arg = blup_rank$Genotyp,
  las = 2,
  col = colors_blup,
  ylab = "Seed Quality (EBLUP)",
  main = "Seed Quality Ranking (BLUP-based)",
  cex.names = 0.8
)


