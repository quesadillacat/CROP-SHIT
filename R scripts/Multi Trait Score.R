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




#ggplot ngl i had a friend help me with this shit

library(ggplot2)
library(dplyr)

highlight_nums <- c(120, 145, 210, 271, 419, 424)

df_plot <- df_ranked %>%
  mutate(
    AssystNum = as.numeric(sub("ASSYST-", "", Genotype)),
    Highlight = ifelse(AssystNum %in% highlight_nums, "Relevant accessions", "Irrelevant accessions")
  )

plot <- ggplot(df_plot, aes(x = Genotype, y = MultiTraitScore, fill = Highlight)) +
  geom_col() +
  scale_fill_manual(values = c(
    "Relevant accessions" = "darkgreen",
    "Irrelevant accessions" = "lightblue"
  )) +
  scale_y_continuous(
    breaks = seq(-6, 6, by = 1),
    limits = c(-6, NA)      
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12, face = "bold"),
    axis.title.x = element_text(size = 14, margin = margin(t = 15)),
    axis.title.y = element_text(size = 14, margin = margin(r = 15))
  ) +
  labs(
    
    x = "Genotype",
    y = "Multi-Trait Score"
  )

ggsave("MultiTraitScore.png", plot = plot, width = 14, height = 7, dpi = 600)

