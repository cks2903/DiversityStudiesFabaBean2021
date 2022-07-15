#install.packages("pcadapt")
library(pcadapt)

setwd("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220228_PCAdapt")

# create input for pcadapt function
###################################

# individuals in columns
# SNPs in lines
# SNPs coded as 0,1,2
path_to_file <- "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220228_PCAdapt/genotypes_1MAFfilter.bed"
filename <- read.pcadapt(path_to_file, type = "bed",type.out="bed") #convert to bed format needed for pcadapt functions
str(filename) # fits with expected number of individuals and SNPs



# Choosing the number K of Principal components
###################################

# set a large K
x <- pcadapt(input = filename, K = 20)

# scree plot to display percentage of var. explained by each PC
plot(x, option = "screeplot")
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. We recommend to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
# K is around 4, I would say


# alternatively we can trust in the faststructure K=4 and define populations and look which PC axises can no longer retain structure
ADMIXTURE = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220224_ADMIXTURE/genotypesADMIXTURE_MAFfiltered.3.Q")
colnames(ADMIXTURE) = paste("Pop",seq(1:ncol(ADMIXTURE)),sep="")
ADMIXTURE$Line= Ind



POP = as.character(rep(0,685))

for (i in seq(1:nrow(ADMIXTURE))){
  belongto = colnames(ADMIXTURE)[which(ADMIXTURE[i,1:ncol(ADMIXTURE)]==max(ADMIXTURE[i,1:ncol(ADMIXTURE)]))]
  if (ADMIXTURE[i,which(ADMIXTURE[i,1:ncol(ADMIXTURE)]==max(ADMIXTURE[i,1:ncol(ADMIXTURE)]))]>0.5){
    POP[i] = belongto
  } else {
   POP[i] = NA
  }
}

poplist.names = POP

plot(x, option = "scores", i= 1, j= 2, pop = poplist.names) # grouping
plot(x, option = "scores", i= 2, j= 3, pop = poplist.names) # grouping
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names) # grouping
plot(x, option = "scores", i = 4, j = 5, pop = poplist.names) # somewhat grouping
plot(x, option = "scores", i = 5, j = 6, pop = poplist.names) # grouping gets less clear
plot(x, option = "scores", i = 6, j = 7, pop = poplist.names) # grouping gets less clear
plot(x, option = "scores", i = 7, j = 8, pop = poplist.names) # grouping gets less clear
plot(x, option = "scores", i = 9, j = 10, pop = poplist.names) # hard to distinguish groups
plot(x, option = "scores", i = 11, j = 12, pop = poplist.names) # hard to distinguish groups
plot(x, option = "scores", i = 13, j = 14, pop = poplist.names) # hard to distinguish groups
plot(x, option = "scores", i = 15, j = 16, pop = poplist.names) # hard to distinguish groups
plot(x, option = "scores", i = 17, j = 18, pop = poplist.names) # hard to distinguish groups
plot(x, option = "scores", i = 19, j = 20, pop = poplist.names) # mixing well
# hard to see could be somewhere from K=4 to K=8 really




# Computing the test statistic based on PCA
###################################
x <- pcadapt(filename, K = 3,min.maf= 0.00)
summary(x)
plot(x , option = "manhattan")
plot(x, option = "qqplot")
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange") # most p-values follow a uniform distribution however we have more small p-values (outliers) and large p-values ( I guess SNPs that do not differ)
plot(x, option = "stat.distribution")

# Choosing a cut of for outlier detection
###################################

# method 1: q-values
#BiocManager::install(c("qvalue"))
library(qvalue)

qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05 #(Like Bayescan, so they are comparable)
outliers <- which(qval < alpha)
length(outliers) # 407

mapfile=read.table("genotypes_1MAFfilter.map",head=F,sep="\t")
outlierSNPs= mapfile[outliers,2]
outlierSNPs

# make table
table=cbind(mapfile[,2],qval,x$maf)
colnames(table)=c("SNP name","qval","MAF")

write.table(table,paste("PCAdaptResults_K3_PopMembershipAbove50",Sys.Date(),".txt",sep="_"),col.names=T,row.names = F,quote=F,sep="\t")

# method 2: Bonferoni correction
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers) # 81



# method 3: Benjamini Hochberg
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers) # 407


snp_pc <- get.pc(x, outliers)


# should account for LD in genome scans, but there is no LD to worry about really
par(mfrow = c(2, 2))
for (i in 1:4){
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
}
