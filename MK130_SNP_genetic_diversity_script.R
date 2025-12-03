###########################################################################
# Part 0
# LOAD packages
###########################################################################

# install package "adegenet" and package "hierfstat"

# set directory - please adapt
setwd("D:/Documents/Genetik_der_Nutzpflanzendiversität/Teaching/MP128_GOCD/WS25_26/")

# load package
library("adegenet")
library("hierfstat")

###########################################################################
# Part 1
# LOAD and FORMAT data
###########################################################################

snp<-read.table(file="MK130_60k_ASSYST_selected.txt", header =T)

# we need to transpose the data for the package
t_snp<-data.frame(t(snp))

# rename
colnames(t_snp)<-t_snp[1,]

# remove the name line
t_snp<-t_snp[-1,]

# extract individual names
ind<-as.character(rownames(t_snp))
# extract population names
population<-as.character(t_snp$type)

# remove the type 
t_snp<-t_snp[,-1]

####to reduce computing time and allow easy display, we only look at Chromosome A01:
t_snp<-t_snp[,-c(2674:52157)]

# convert to genind object
snp_data<-df2genind(t_snp, ploidy=2, ind.names=ind, pop=population, sep="")
# look at it
snp_data

###########################################################################
# Part 2
# CALCULATE  and DISPLAY genetic distance between accessions
###########################################################################
# Create hierfstat object first

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


###########################################################################
# Part 3
# CALCULATE  and DISPLAY measures of genetic diversity
###########################################################################

# calculate number of alleles per locus

number_of_alleles<-nAll(snp_data)
mean(number_of_alleles)

plot(number_of_alleles, xlab="Loci number", ylab="Number of alleles", 
     main="Number of alleles per locus")

# calculate allele richness
allrich<-allelic.richness(snp_data,min.n=NULL,diploid=TRUE)
head(allrich$Ar, 5)

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


