rm(list=ls())
set.seed(1234567)
setwd("H:/Final Thesis/Gnome Analysis/R Codes")
source("Functions4GWAS.R")
##########################
## READ DATA ####
##########################
# Read the new files
pop <- read.csv("POP_DNA.csv", sep = ",")
pt <- read.csv("PT.csv", sep = ",")
gene.location <- read.csv("SiteFile.csv", sep = ",")
org_zero <- read.csv("KO_file_0.csv", sep = ",")
org_one <- read.csv("KO_file_1.csv", sep = ",")
gen <- read.csv("GEN.csv", sep = ",")
#lod <- read.csv("LOD.csv", sep = ",")


# Define the size of the data 
N = 999 #Number of individuals
genome.size = 160000 #Length of the genome

##########################################
# Create a genome matrix G with letters
#Gtext <- matrix(0, N, genome.size)
#for (i in 1:N)
#{
#  Gtext[i,] <- unlist(strsplit(substr(pop$genome[i], 1,genome.size), ""))
#}
#pop$genome[1]
#save(Gtext, file = "Gtext.RData")

load("Gtext.RData")

##########################
## CREATE SNP DATA ####
##########################
#SNP_data <- G_to_binary(Gtext) #Recode the matrix Gtext to binary
#save(SNP_data, file = "SNP_data.RData")

load("SNP_data.RData")

###################################################################
## LOOK AT HOW SIMILAR THE GENOMES ARE FOR DIFFERENT INDIVIDUALS####
######################################################################
test1 <- colMeans(SNP_data)<0.95 #MAF > 5 percent
subset1 <- seq(1, sum(test1), length.out=5000) #Enough with 5000 markers to compute relatedness
GR <- Relatedness(SNP_data[,test1][,subset1], PLOT=TRUE)


#library(RColorBrewer)
#heatmap(GR)


##########################
## GET THE PHENOTYPES  ####
##########################


Phenotype <- get.phenotype(pt, trait="scoreIndividual")
boxplot(Phenotype, main="Boxplot for scoreIndividual Phenotype")
hist(Phenotype)


# ***** Check Relationship of features (with Mean) *****

phenotype_pairs <- Phenotype_Mean()

pairs(phenotype_pairs)


# ***** Check Relationship of features (with Variance) *****

phenotype_pairs_var <- Phenotype_Var()

pairs(phenotype_pairs_var)




#####################################
## COMPUTE THE HERITABILITY  ########
#####################################

h2 <- compute.heritability(Phenotype, GR)
cat("Estimated heritability:", h2, "\n")


#########################################################################################
## HOW WELL CAN WE PREDICT INDIVIDUAL PHENOTYPES BASED ON THE GENETIC DATA ONLY ########
## This model assumes all SNP have a small and normally distributed effect on the trait 
## and that these small effects are additive.
#########################################################################################
CP <- cross.predict(Phenotype, GR)
par(mfrow=c(1,1))
plot(y=CP$predicted, x=CP$unobserved.phenotypes,
     main="Predicting Phenotypes",
     xlab="Observed", ylab="Predicted",
     ylim=range(c(CP$predicted,CP$unobserved.phenotypes))
     )
abline(0,1, col="grey")
cat("Correlation between predicted and observed", cor(CP$predicted, CP$unobserved.phenotypes), "\n")



# Creating a Genome matrix G
G <- matrix(0, N, genome.size)
for (i in 1:N)
{
  G[i,] <- unlist(strsplit(substr(pop$genome[i], 1,genome.size), ""))
}

#A function to change letters to 0 or 1
G_to_binary <- function(G) {
  M <- matrix(0, nrow(G), ncol(G))
  for (j in 1:ncol(G)) {
    most.common <- names(which.max(table(G[,j])))
    M[G[,j]==most.common,j] = 1
  }
  return(M)
}

summary(h2)
##########################
## ANALYSIS OF DATA ####
##########################
#summary(Phenotype)
hist(Phenotype)
SNP_data <- G_to_binary(G) #Recode the matrix G to be able to compute correlations
#View(SNP_data)
GWAS_results <- numeric(genome.size)
for (j in 1:genome.size) {
  if (sum(SNP_data[,j]) < (N-5)) {
    mod1 <- lm(Phenotype~SNP_data[,j])
    P_value <- summary(mod1)$coef[2,4]
  } else {
    P_value=1
  }
  GWAS_results[j] <- -log10(P_value)
}


table(G[,which.max(GWAS_results)])
hist(GWAS_results)
plot(GWAS_results, xlab="Nucleotide position", ylab="-log10(P-value)")
median(GWAS_results)


which.max(GWAS_results)
#GWAS_results[28418]

"genomic inflation factor"
# genome-wide association study

indx = (which.max(GWAS_results) - 100): (which.max(GWAS_results) +100)
subset_SNPdata= SNP_data[,indx]


Sub_SNPdata <- GWAS_results[]

#abline(h = 7.3)


# ********** New Analysis ***********

SNP_data = subset_SNPdata
MAF <- 1-colSums(SNP_data)/nrow(SNP_data) #Compute minor allele frequency
plot(MAF)
hist(MAF)
test <- MAF>0
print(sum(test)) #How many positions have varying nucleotide information?
SNP_data <- SNP_data[,test] 
M <- matrix(NA, ncol(SNP_data), ncol(SNP_data))
#mm = matrix(NA, ncol(SNP_data)^2, 3)
k=0
for ( i in 1:(ncol(SNP_data)-1) ) {
  for (j in (i+1):ncol(SNP_data)) {
    k=k+1
    M[i,j] <- cor(SNP_data[,i], SNP_data[,j])
    
  }
}

M1 <- M
M1[is.na(M)]=0
M2<- M1+t(M1)
diag(M2)=1
n <- ncol(M2)
step = 100
indx = step
kk <- ceiling(sqrt(n/step))
diag(M2)=NA
library(RColorBrewer)
pdf()
while(indx < n) {
  heatmap(M2[(indx-step+1):indx,(indx-step+1):indx], Colv=NA, Rowv=NA)
  legend(x="bottomright", legend=c("min", "ave", "max"), 
         fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))
  indx = indx +step
}
heatmap(M2[(indx-step+1):n,(indx-step+1):n], Colv=NA, Rowv=NA)
legend(x="bottomright", legend=c("min", "ave", "max"), 
       fill=colorRampPalette(brewer.pal(8, "Oranges"))(3))
dev.off()


# ********************* 

lm(Phenotype~SNP_data[,1])
mod1 = lm(Phenotype~SNP_data[,1])
summary(mod1)
cor.test(Phenotype, SNP_data[,1])