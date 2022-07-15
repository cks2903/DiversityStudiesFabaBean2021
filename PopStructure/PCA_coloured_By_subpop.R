library(tidyverse)
library(viridis)
library(plyr)

setwd("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220224_ADMIXTURE")
pca <- read_table2("pca.eigenvec", col_names = FALSE)
eigenval <- scan("pca.eigenval")

# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
pve <- data.frame(PC = 1:nrow(pca), pve = eigenval/sum(eigenval)*100)

# get panel information
Subpop3<- read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220224_ADMIXTURE/genotypesADMIXTURE_MAFfiltered.3.Q",head=F,sep="\t")
colnames(Subpop3)=c("Subpop1","Subpop2","Subpop3")
indnames = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220224_ADMIXTURE/BedFileForStructureFiltered_20220223.fam",head=F,sep=" ")
SubpopInfo = cbind(indnames[,1],Subpop3)
colnames(SubpopInfo)[1]="Accession"

SubpopInfo$POP = as.character(rep(0,nrow(SubpopInfo)))

for (i in seq(1:nrow(SubpopInfo))){
  belongto = colnames(SubpopInfo)[which(SubpopInfo[i,1:ncol(SubpopInfo)]==max(SubpopInfo[i,2:(ncol(SubpopInfo)-1)]))]
  membershipcoef = SubpopInfo[i,which(colnames(SubpopInfo)==belongto)]
  if (membershipcoef>=0.5){
    SubpopInfo$POP[i] = belongto
  } else {
    SubpopInfo$POP[i] = NA
  }
}

# now determine which pop it a line belongs to


pca_copy = as.data.frame(pca)
colnames(SubpopInfo)[1]="ind"
merged <- merge(pca_copy,SubpopInfo,by="ind",sort = FALSE)
all(merged[,1]==pca[,1])

# include geographic origin
Origin <- read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220224_ADMIXTURE/Origin_Broad.txt",sep="\t",head=T)
colnames(Origin)[1]="ind"
SubpopInfo_merged = join(SubpopInfo,Origin)
all(SubpopInfo_merged[,1]==SubpopInfo[,1]) #check


#remake
pca <- as.tibble(data.frame(pca, as.character(SubpopInfo$POP),SubpopInfo_merged$OriginLessBroad))


#plotting


# plot pca
min= min(c(min(pca$PC1),min(pca$PC2)))
max =max(c(max(pca$PC1),max(pca$PC2)))


c <- ggplot(pca, aes(PC1, PC2, col =SubpopInfo.POP)) + geom_point(size = 3)
#c <- c +   scale_colour_manual(values = c("#0D2C54","#00AF54","#C6878F")) 
c <- c + coord_equal() + theme_classic() + xlim(min,max) +ylim(min,max) +scale_color_manual(values = c("#440154","#fde725","#21918c")) 
d<-c + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
ggsave("PC1PC2_3subpops.png", plot = d, width = 20, height = 15, unit = 'cm')


pca$SubpopInfo_merged.OriginLessBroad[is.na(pca$SubpopInfo_merged.OriginLessBroad)]="No origin"
c <- ggplot(pca, aes(PC1, PC2, col =SubpopInfo.POP,shape = SubpopInfo_merged.OriginLessBroad)) + geom_point(size = 3)
#c <- c +   scale_colour_manual(values = c("#0D2C54","#00AF54","#C6878F")) 
c <- c + coord_equal() + theme_classic() + xlim(min,max) +ylim(min,max) +scale_color_manual(values = c("#440154","#fde725","#21918c")) +  scale_shape_manual(values = c(0, 2,5,4,1))
d<-c + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
ggsave("PC1PC2_3subpops_colouredByOrigin.png", plot = d, width = 20, height = 15, unit = 'cm')



e <- ggplot(merged, aes(PC2, PC3, col =POP )) + geom_point(size = 1)
#e <- e +   scale_colour_manual(values = c("#0D2C54","#00AF54","#C6878F")) 
e <- e + coord_equal() + theme_classic()  + xlim(min,max) +ylim(min,max) +scale_color_manual(values = c("#440154","#fde725","#21918c")) 
f<-e + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))
ggsave("PC2PC3_3subpops.png", plot = f, width = 20, height = 15, unit = 'cm')


c <- ggplot(pca, aes(PC2, PC3, col =SubpopInfo.POP,shape = SubpopInfo_merged.OriginLessBroad)) + geom_point(size = 3)
#c <- c +   scale_colour_manual(values = c("#0D2C54","#00AF54","#C6878F")) 
c <- c + coord_equal() + theme_classic() + xlim(min,max) +ylim(min,max) +scale_color_manual(values = c("#440154","#fde725","#21918c")) +  scale_shape_manual(values = c(0, 2,5,4,1))
d<-c + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))
ggsave("PC2PC3_3subpops_colouredByOrigin.png", plot = d, width = 20, height = 15, unit = 'cm')



