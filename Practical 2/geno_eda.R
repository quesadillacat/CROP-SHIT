library(readxl)
library(dplyr)
library(tidyr)
library(adegenet)
library(hierfstat)

kasp_dec3_raw <- read_excel("KASP shit/2025-12-03_081443.xlsx")
kasp_dec4_raw <- read_excel("KASP shit/2025-12-04_084604.xlsx")
kasp_mk130_dna <- read_excel("KASP shit/MK130_DNA_Nov2025_second_expierments_011_27_2025.xlsx")

pflzn <- data.frame(
  accessions = c("BnASSYST_120", "BnASSYST_145", "BnASSYST_210", "BnASSYST_271", "BnASSYST_419", "BnASSYST_424"),
  name       = c("Darmor",       "Hokkai 3-Go",  "RED RUSSIAN",  "Liho",         "Angus",        "Conqueror Bronze Green Top"),
  type       = c("Winter OSR",   "Winter OSR",   "Siber. Kale",  "Spring OSR",   "swede",        "swede")
  )

# Read in SNP chip data, transpose, give column names from 1st row and trim
snp_raw <- read.table("MK130_60k_ASSYST_selected.txt", header = T)
snp <- data.frame(t(snp_raw))
colnames(snp) <- snp[1,]
snp <- snp[-1,]


# Set accessions as row IDs, extract those as our individuals
# Extract crop types as populations
# indivs <- as.character(rownames(snp))
# populations <- as.character(snp$type)



# to only grab our accessions
snp_trim <- snp %>%
  filter(row.names(snp) %in% pflzn$accessions)

indivs <- as.character(rownames(snp_trim))
populations <- as.character(snp$type)
pflzn %>%
  mutate(pflzn$type = populations)

# trim off type (we'll use these extracted fields for genind instantiation)
snp_2_genind <- snp_trim %>%
  select(!type)











###############################################################################
# tidyverse re-write for my own readability reasons
###############################################################################

snp.ours <- snp_raw %>%
  select(pflzn$accessions)

indivs.ours <- pflzn$accessions
populations.ours <- snp.ours