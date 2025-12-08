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
populations <- as.character(snp_trim$type)
pflzn %>%
  mutate(pflzn$type = populations)

# trim off type (we'll use these extracted fields for genind instantiation)
snp_2_genind <- snp_trim %>%
  select(!type)

#### to reduce computing time and allow easy display, we only look at Chromosome A01:
# snp_a01 <- snp_2_genind[,-c(2674:52157)]
snp_wholegenome <- snp_2_genind
rm(snp_2_genind, snp_trim) # optional: rm intermediate data objects we're not using anymore

# convert to genind
snp_data<-df2genind(snp_wholegenome, ploidy=2, ind.names=indivs, pop=populations, sep="")


# GENETIC DISTANCE PCA PLOT
snp_data_1 <- genind2hierfstat(snp_data)

# calculate PCA
x<-indpca(snp_data_1, ind.labels = rownames(snp_data_1))

popcol<-c(rep("blue",3), rep("turquoise",3), rep("black",3),rep("purple",3))
# plot PCA1/2
plot(x, cex=0.8, col=popcol, ax1=1, ax2=2)
# plot PCA1/3
plot(x, cex=0.8, col=popcol,ax1=1, ax2=3)
# plot PCA2/3
plot(x, cex=0.8, col=popcol,ax1=2, ax2=3)

# CALC GENETIC DIVERSITY

# No. alleles per locus

allele.num <- nAll(snp_data)
allele.num.avg <- mean(number_of_alleles)

plot(number_of_alleles, xlab="Loci number", ylab="Number of alleles", 
     main="Number of alleles per locus")

# calculate allele richness
allele.rich <- allelic.richness(snp_data,min.n=NULL,diploid=TRUE)
head(allrich$Ar, 5) # head() call unnecessary?

# calculate observed and expected heterozygosity
# using adegenet
div<-summary(snp_data)
names(div)

plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")

plot(div$Hexp, xlab="Loci number", ylab="Expected Heterozygosity", 
     main="Expected heterozygosity per locus")

# calculate basic statistics from hierfstat
basicstat <- basic.stats(snp_data_1, diploid = TRUE, digits = 2)

# extract locus-specific data
locus_results<-basicstat$perloc
head(locus_results)

# plot observed heterozygosity 
# (and compare against adegenet results in blue)
plot(locus_results$Ho, xlab="Loci number", ylab="Ho", 
     main="Observed heterozygosity per locus")
points(div$Hobs, col="blue")

# plot expected heterozygosity within subpopulations
plot(locus_results$Hs, xlab="Loci number", ylab="Hs", 
     main="Expected heterozygosity per locus within subpopulations")

# plot expected heterozygosity in the total population 
# (and compare against adegenet results in blue)
plot(locus_results$Ht, xlab="Loci number", ylab="Ht", 
     main="Expected heterozygosity per locus")
points(div$Hexp, col="blue")

#plot Dst (Ht-Hs) per locus
plot((locus_results$Dst), xlab="Loci number", ylab="Dst", 
     main="Subpopulation differentiation")

#plot Dstp (Ht-Hs corrected) per locus
plot((locus_results$Dstp), xlab="Loci number", ylab="Dst", 
     main="Subpopulation differentiation, corrected")

# plot Fst per locus
plot(locus_results$Fst, xlab="Loci number", ylab="Fst", 
     main="Fst per locus")

# plot Fstp per locus
plot(locus_results$Fstp, xlab="Loci number", ylab="Fst", 
     main="Fst corrected per locus")

# plot Fis per locus
plot(locus_results$Fis, xlab="Loci number", ylab="Fis", 
     main="Fis per locus")

# plot Jost's D per locus
plot(locus_results$Dest, xlab="Loci number", ylab="Jost's D", 
     main="Jost's D per locus")


# compare results corrected/uncorrected
plot((locus_results$Dst), xlab="Loci number", ylab="Dst", 
     main="Subpopulation differentiation", ylim=c(-1,1))
points(locus_results$Dstp, col="blue")


plot(locus_results$Fst, xlab="Loci number", ylab="Fst", 
     main="Fst per locus")
points(locus_results$Fstp, col="green")












###############################################################################
# tidyverse re-write for my own readability reasons
###############################################################################

snp.ours <- snp_raw %>%
  select(pflzn$accessions)

indivs.ours <- pflzn$accessions
populations.ours <- snp.ours