##################################################################################
############ Expected heterozygosity                                  ############ 
##################################################################################

# load libraries needed
#install.packages("inbreedR")
library("adegenet")
library(data.table)

setwd("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20221106_Heterozygosity/")


# load functions
{
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
}

# load data
{
  genotypes=fread("genotypes_filtered.recode.vcf") #vcf file, not MAF filtered
  dim(genotypes)
  genotypes = as.data.frame(genotypes)
  
  
  # Pop info
  SP1inds = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20221106_Heterozygosity/population_1_50perc.txt")[,1]
  SP2inds = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20221106_Heterozygosity/population_2_50perc.txt")[,1]
  SP3inds = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20221106_Heterozygosity/population_3_50perc.txt")[,1]
  
  individuals = which(colnames(genotypes) %in% SP1inds | colnames(genotypes) %in% SP2inds  | colnames(genotypes) %in% SP3inds)
  genotypes_shorter=genotypes[,individuals]
  genotypes_t=t(genotypes_shorter)
  dim(genotypes_t)
  rownames(genotypes_t)=colnames(genotypes_shorter)
  colnames(genotypes_t)=genotypes[,3]
  genotypes=genotypes_t
  names=rownames(genotypes)
  genotypes[genotypes=="1|1"] = "1 1"
  genotypes[genotypes=="1|0"] = "0 1"
  genotypes[genotypes=="0|1"] = "0 1"
  genotypes[genotypes=="0|0"] = "0 0"
}


Pop=c()
for (i in seq(1:nrow(genotypes))){
  
  currentind=rownames(genotypes)[i]
  
  if (currentind %in% SP1inds){
    
    Pop=c(Pop,"SP1")
  }
  if (currentind %in% SP2inds){
    
    Pop=c(Pop,"SP2")
  }
  if (currentind %in% SP3inds){
    Pop=c(Pop,"SP3")
    
  }
}


POP = as.factor(Pop)

obj_geno <- df2genind(genotypes, ploidy=2, sep=" ")
obj_geno
obj_geno@pop=POP



# calculate heterozygosity
summ= summary(obj_geno)
Exp_het=Hs(obj_geno) # expected heterozygosity within populations

# Observed per pop
n.pop  <- seppop(obj_geno) 
mean.hobs <- do.call("c", lapply(n.pop, function(x) mean(summary(x)$Hobs))) 


# calculate observed heterozygosity in one panel for instance
write.table(mean.hobs, "ObservedHeterozygosity.txt",sep="\t",quote=F)

write.table(Exp_het, "ExpectedHeterozygosity.txt",sep="\t",quote=F)
