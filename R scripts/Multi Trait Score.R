# Multi trait score!!

# Load dataset containing one predicted value (EBLUP) per genotype and trait
df <- read.table("MTS.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE)

# Convert all trait columns to numeric
trait_columns <- c("SeedQuality", "CCM", "HGW", "DW", "FM", "PH", "BBCH")
df[trait_columns] <- lapply(df[trait_columns], as.numeric)


# Standardize traits so that all traits are on the same scale

# Standardization  ensures that no single trait dominates
# the multi-trait score just because it has larger numerical values so this is totally critical
# Seed Quality was already standardized earlier when I calculated Seed quality index earlier

traits_to_scale <- c("CCM", "HGW", "DW", "FM", "PH", "BBCH")
df[traits_to_scale] <- scale(df[traits_to_scale])



# The multi-trait score is the sum of all standardized trait values
# for each genotype.

# High score= strong performance across many traits
# Low score  = weaker overall performance.

df$MultiTraitScore <- rowSums(df[, c("SeedQuality", traits_to_scale)]) # adds everything together


# Rank genotypes by their multi-trait score high to low

df_ranked <- df[order(-df$MultiTraitScore), ]

# Extract accession numbers
# and sort by number

df_ranked$AssystNum <- as.numeric(sub("ASSYST-", "", df_ranked$Genotype))
df_sorted_by_number <- df_ranked[order(df_ranked$AssystNum), ]

# aaand extract winners for each trait
traits <- c("SeedQuality", "CCM", "HGW", "DW", "FM", "PH", "BBCH")
winners <- sapply(
  traits,
  function(trait_name) {
    df_ranked$Genotype[ which.max(df_ranked[[trait_name]]) ]
  }
)


df_ranked
df_sorted_by_number
winners


