##################################################################################
############ Expected heterozygosity                                  ############ 
##################################################################################

# load libraries needed
#install.packages("inbreedR")
library("adegenet")

args=commandArgs(trailingOnly = TRUE)


# load functions
{
  is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
}

# load data
{
genotypes=read.table("/home/cks/norfab/faststorage/cks/20220120_DiversityPanel/20220120_Panel_Specific_Analyses/20220126_RedundancyWithinPanels/Mapped_HighQualSNPs2022_RedFiltered20220126.vcf", header = T, sep = "\t",stringsAsFactors = F) #vcf file, not MAF filtered
dim(genotypes)
genotypes_shorter=genotypes[,10:ncol(genotypes)]
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


# Panel info
Dic <- read.table("/home/cks/norfab/faststorage/cks/20220120_DiversityPanel/20220120_Panel_Specific_Analyses/20220126_RedundancyWithinPanels/Lines_Panel_DicNew.txt",head=T,sep="\t")
Panels = unique(Dic[,2])
X = Dic[which(Dic[,1] %in% rownames(genotypes)),1:2]

# Sort X after genotypes
order(X[,1])

X=X[order(match(X[,1],rownames(genotypes))),]
all(X[,1]==rownames(genotypes))

X[which(X[,2] %in% c("VICCI_O1","VICCI_O4_Sel","VICCI_O2_Un","VICCI_O6_Un","VICCI_O4_Un","VICCI_O2_Sel","VICCI_O6_Sel")),2]="VICCI"


POP = as.factor(X[,2])

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
write.table(mean.hobs, "/home/cks/norfab/faststorage/cks/20220120_DiversityPanel/20220120_Panel_Specific_Analyses/20220126_Heterozygosity/ObservedHeterozygosity.txt",sep="\t",quote=F)

write.table(Exp_het, "/home/cks/norfab/faststorage/cks/20220120_DiversityPanel/20220120_Panel_Specific_Analyses/20220126_Heterozygosity/ExpectedHeterozygosity.txt",sep="\t",quote=F)
