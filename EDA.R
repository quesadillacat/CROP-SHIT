library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

pheno <- read_excel("Pheno Data.xlsx")

df <- pheno %>%
  group_by(Genotype) %>%
  mutate(Rep = row_number()) %>%
  ungroup()