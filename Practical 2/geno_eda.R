library(readxl)
library(dplyr)


kasp_dec3_raw <- read_excel("KASP shit/2025-12-03_081443.xlsx")
kasp_dec4_raw <- read_excel("KASP shit/2025-12-04_084604.xlsx")
kasp_mk130_dna <- read_excel("KASP shit/MK130_DNA_Nov2025_second_expierments_011_27_2025.xlsx")

our.plants <- c("ASSYST-120", "ASSYST-145", "ASSYST-210", "ASSYST-271", "ASSYST-419", "ASSYST-424")

pflzn <- data.frame(
  accessions = c("ASSYST-120", "ASSYST-145", "ASSYST-210", "ASSYST-271", "ASSYST-419", "ASSYST-424"),
  name       = c("Darmor"),
  type       = c("WOSR")
  , row.names = FALSE, stringsAsFactors = TRUE)