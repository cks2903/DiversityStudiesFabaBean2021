
#### Load libraries ####
library("pegas")
library('adegenet')
library("poppr") #perhaps use?




genotype=read.table("genotypes_1MAFfilter.ped",head=F,sep=" ",stringsAsFactors=F) # ped file from PLINK
dim(genotype)
genotype_only=genotype[,7:ncol(genotype)]

SNPS=read.table("genotypes_1MAFfilter.map",head=F,sep="\t") # map file from PLINK
SNPnames= SNPS[,2]

structure=read.csv("ForMapping.csv",sep=",",head=T) # a file that contains samples in first col and SP membership in second column
head(structure)


# make sure it is in same order as geno1, and remove the 20 individuals not in any pop. from genotype_only
famtile=read.table("genotypes_1MAFfilter.fam",head=F,sep=" ") # fam file from PLINK
dim(famtile)

order_of_ped_files= famtile[,1]
all(order_of_ped_files==structure$Sample) #order is not the same

structure_ordered=structure[ order(match(structure$Sample, order_of_ped_files)), ]
all(order_of_ped_files==structure_ordered$Sample) #order is now the same

population_hier=structure_ordered[,c(1,ncol(structure_ordered))]
colnames(population_hier)=c("Accession","Pop")
rownames(population_hier)=population_hier$Accession
population_hier$Accession=NULL

remove=which(is.na(structure_ordered$Pop)=="TRUE")
structure_ordered=structure_ordered[-remove,]
genotype_only=genotype_only[-remove,]


# convert to the required format (df2genind function works on  dataframe w. genotypes are one row pr. genotype, markers in columns is element is a string of coding  characters)
geno1=df2genind(genotype_only,sep=" ")
geno1 # genid format

#define your population groups
geno1@pop = as.factor(population_hier$Pop)
geno1@other$Pop= as.factor(population_hier$Pop)
g=as.factor(structure_ordered$Pop)

# Do amova
geno1_dist  <- dist(geno1) # Euclidean distance
amova <- pegas::amova(geno1_dist ~g,nperm = 10000) # 1 level
amova
str(amova)

#save results from AMOVA
save(amova, file=paste("amova",Sys.Date(),".Rdata"))
load("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220314_AMOVA/amova 2022-03-15 .Rdata")


# poppr
geno2=geno1
geno2@pop=g
other(geno2)$population_hierarchy=g
strata(geno2) <-  data.frame(other(geno2)$pop)


agc <- as.genclone(geno2)
agc

amova.result <- poppr.amova(agc, ~other.geno2..pop)

