#####################################################################
#       make BayeScan SNP input file, with K=3 and pop. membership is 50%
#####################################################################

library("varhandle")

# 
mapfile = read.table("MAF1perc_filtered.map",head=F, sep="\t")

FarmCPU_genotypeinfo =read.table("/home/cks/norfab/faststorage/cks/20220120_DiversityPanel/20220126_DiversityPanel/RelationshipMatrix/GenotypeInfo_FarmCPU.txt",head=T, sep="\t") 
FarmCPU_genotypes =read.table("/home/cks/norfab/faststorage/cks/20220120_DiversityPanel/20220126_DiversityPanel/RelationshipMatrix/Genotypes_FarmCPU.txt",head=T, sep="\t",check.names = F) 

# filter away SNPs with too low MAF <0.01, that is not in map file
keep <- which(FarmCPU_genotypeinfo$Name %in% mapfile[,2])
all(FarmCPU_genotypeinfo[,1] == colnames(FarmCPU_genotypes)[2:ncol(FarmCPU_genotypes)]) # be sure

FarmCPU_genotypeinfo = FarmCPU_genotypeinfo[keep,]
keep2 <- which(colnames(FarmCPU_genotypes) %in% FarmCPU_genotypeinfo$Name)
FarmCPU_genotypes = FarmCPU_genotypes[,c(1,keep2)]



# one row for each ind
nrow(FarmCPU_genotypes) #fits

# individual names
Ind=FarmCPU_genotypes[,1]



# merge with pop information on each individual

ADMIXURE = read.table("/home/cks/norfab/faststorage/cks/20220120_DiversityPanel/20220126_DiversityPanel/STRUCTURE/genotypesADMIXTURE_MAFfiltered.3.Q")
#Faststructure_d = read.table("structure.6.meanQ")
colnames(ADMIXURE) = paste("Pop",seq(1:ncol(ADMIXURE)),sep="")
ADMIXURE$Line= Ind

FarmCPU_genotypes$POP = as.character(rep(0,nrow(FarmCPU_genotypes)))

for (i in seq(1:nrow(ADMIXURE))){
  belongto = colnames(ADMIXURE)[which(ADMIXURE[i,1:(ncol(ADMIXURE)-1)]==max(ADMIXURE[i,1:(ncol(ADMIXURE)-1)]))]
  if (ADMIXURE[i,which(ADMIXURE[i,1:(ncol(ADMIXURE)-1)]==max(ADMIXURE[i,1:(ncol(ADMIXURE)-1)]))]>0.5){
    FarmCPU_genotypes$POP[i] = belongto
  } else {
    FarmCPU_genotypes$POP[i] = NA
  }
}

# rearrange to desired format
Desired = FarmCPU_genotypes[,c(21118,1:21117)]
DesiredFormat =cbind(Ind,Desired)

# now replace 
DesiredFormat$POP[which(DesiredFormat$POP=="Pop1")]=as.character(1)
DesiredFormat$POP[which(DesiredFormat$POP=="Pop2")]=as.character(2)
DesiredFormat$POP[which(DesiredFormat$POP=="Pop3")]=as.character(3)

# THEN sort after population
DesiredFormat_sorted <- DesiredFormat[order(DesiredFormat$POP),]
DesiredFormat_sorted$Ind
POP1idx=which(DesiredFormat_sorted$POP=="1")
POP2idx=which(DesiredFormat_sorted$POP=="2")
POP3idx=which(DesiredFormat_sorted$POP=="3")

NoPOPidx=which(is.na(DesiredFormat_sorted$POP)==TRUE)

DesiredFormat_sorted$IndNum=unfactor(DesiredFormat_sorted$Ind)

DesiredFormat_sorted$IndNum[POP1idx]=as.character(seq(1:length(POP1idx)))
DesiredFormat_sorted$IndNum[POP2idx]=as.character(seq(1:length(POP2idx)))
DesiredFormat_sorted$IndNum[POP3idx]=as.character(seq(1:length(POP3idx)))
DesiredFormat_sorted=DesiredFormat_sorted[-NoPOPidx,]

DesiredFormat_sorted1=DesiredFormat_sorted[,c(21120,seq(1:21119))]


# then reset individual name in each pop
DesiredFormat_sorted_2=DesiredFormat_sorted1[,c(1,3:21120)]
DesiredFormat_sorted_2 = DesiredFormat_sorted_2[,-3]
# write file
write.table(DesiredFormat_sorted_2,"BayeScan_input_K=3_Membership50Perc.txt",quote=F,col.names=F,row.names = F)

write.table(DesiredFormat_sorted1,"Individual_numbers_K=3_Membership50Perc.txt",quote=F,col.names=F,row.names = F)
