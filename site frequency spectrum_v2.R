library(ggplot2)

#################################################
#################################################
##    Calculate counts of non-reference alleles
#################################################
#################################################

# load file, not imputed, not filtered for MAF
vcf <- read.table("../Mapped_HighQualSNPs2021-10-01.vcf",sep="\t",header=T)

# List of individuals and the panel they belong to
origin=read.table("/home/cks/norfab/faststorage/cks/20211004_DiversityPanel/20211004_Panel_Specific_Analyses/20211004_Heterozygosity/Lines_Panel_Dic.txt",header=T,sep="\t",na.strings="\\N", quote="",stringsAsFactors = F)
head(origin)
dim(origin)

# for one panel at a time
origin$Panel.ID[which(origin$Panel.ID %in% c("VICCI_O1","VICCI_O4_Sel","VICCI_O2_Un","VICCI_O6_Un","VICCI_O4_Un","VICCI_O2_Sel","VICCI_O6_Sel"))]="VICCI"

Panels <- unique(origin$Panel.ID)

for (panel in Panels){
  
  Accessions_to_keep <- origin$probeset_id[which(origin$Panel.ID==panel)]
  Idx_to_keep <- which(colnames(vcf) %in% Accessions_to_keep)
  vcf_temp <- vcf[,c(seq(1:9),Idx_to_keep)]
  vcf_short <-vcf_temp[,10:ncol(vcf_temp)]
  
  # now start counting how often the reference allele occur for each SNP
  AltAlleleCount = data.frame(matrix(NA, nrow = nrow(vcf_temp), ncol = 3))
  colnames(AltAlleleCount) = c("SNP","Ref_counts","Alt_counts")
  for (i in seq(1:nrow(vcf_temp))){
    
    alt_allele = vcf_temp$ALT[i]
    ref_allele = vcf_temp$REF[i]
    
    homozygote_alt = paste(alt_allele,"/",alt_allele,sep="")
    homozygote_ref = paste(ref_allele,"/",ref_allele,sep="")
    heterozygotes1 =   paste(alt_allele,"/",ref_allele,sep="")
    heterozygotes2 =   paste(ref_allele,"/",alt_allele,sep="")
    
    homozygote_alt_count <- length(which(vcf_short[i,]==homozygote_alt))
    homozygote_ref_count <- length(which(vcf_short[i,]==homozygote_ref))
    heterozygote_count <- length(which(vcf_short[i,]==heterozygotes1 | vcf_short[i,]==heterozygotes2))
    Missing <- length(which(vcf_short[i,]=="./."))
    
    #check
    check1 <- homozygote_alt_count + homozygote_ref_count+ heterozygote_count + Missing == length(Accessions_to_keep)
    if (check1 != TRUE){
      print("There is a problem at check1")
    }
    
    Alt_allelefreq <- round((homozygote_alt_count+heterozygote_count*0.5)/length(Accessions_to_keep),2)
    ref_allelefreq <- round((homozygote_ref_count+heterozygote_count*0.5)/length(Accessions_to_keep),2)
    missingfreq <- round(Missing/length(Accessions_to_keep),2)
    
    check2 <- Alt_allelefreq+ref_allelefreq+missingfreq
    
    if (check2 < 0.989){
      print(paste("There is a problem at check2",i))
    }
    AltAlleleCount$SNP[i] = as.character(vcf_temp$ID[i])
    AltAlleleCount$Ref_counts[i] = homozygote_ref_count*2 + heterozygote_count
    AltAlleleCount$Alt_counts[i] = homozygote_alt_count*2 + heterozygote_count

  }
  
  
  write.table(AltAlleleCount,paste("AlternativeAlleleCount",panel,".txt"),col.names = T,row.names = F,quote=F,sep="\t")
  
  # plotting with non-polymorphic SNPs
  p1<-ggplot(AltAlleleCount, aes(x=Alt_counts)) + 
    geom_histogram(aes(y = ..count..),fill="black",binwidth=1,alpha=0.7)+
    xlab("Alternative allele count") +
    ylab("Number of SNPs") +
    theme_classic()
  p1
  filename_p1= paste("SFS_AllSNPs",panel,Sys.Date(),".pdf",sep="")
  ggsave(filename_p1, plot = p1, width = 20, height = 14, unit = 'cm')
  
  # plotting only polymorphic SNPs
  rmv = which(AltAlleleCount$Alt_counts==0 | AltAlleleCount$Alt_counts==ncol(vcf_short)*2)
  AltAlleleCount_filt = AltAlleleCount[-rmv,]
  
  p15<-ggplot(AltAlleleCount_filt, aes(x=Alt_counts)) + 
    geom_histogram(aes(y = ..count..),fill="black",binwidth=1,alpha=0.7)+
    xlab("Alternative allele count") +
    ylab("Number of SNPs") +
    theme_classic()
  p15
  filename_p15= paste("SFS_AllPolymorphicSNPs_",panel,Sys.Date(),".pdf",sep="")
  ggsave(filename_p15, plot = p15, width = 20, height = 14, unit = 'cm')
  

  # folded SFS
  AltAlleleCount_filt$MAC = AltAlleleCount_filt$Alt_counts
  for (i in seq(1:nrow(AltAlleleCount_filt))){
    if (AltAlleleCount_filt$Alt_counts[i]> length(Accessions_to_keep)){
      AltAlleleCount_filt$MAC[i] =  AltAlleleCount_filt$Alt_counts[i] -length(Accessions_to_keep)
    } 
  }
  p2<-ggplot(AltAlleleCount_filt, aes(x=MAC)) + 
    geom_histogram(aes(y = ..count..),fill="black",binwidth=1,alpha=0.7)+
    xlab("Alternative allele count") +
    ylab("Number of SNPs") +
    theme_classic()
  p2
  
  filename_p2= paste("FoldedSFS",panel,Sys.Date(),".pdf",sep="")
  ggsave(filename_p2, plot = p2, width = 20, height = 14, unit = 'cm')
}






