library(ggplot2)
library(varhandle)
args=commandArgs(trailingOnly = TRUE)

#################################################
#################################################
##    Calculate counts of non-reference alleles
#################################################
#################################################

# load file, not imputed, not filtered for MAF
vcf <- read.table(args[1],sep="\t",header=F)

vcf_short <-vcf[,10:ncol(vcf)]
  
# now start counting how often the reference allele occur for each SNP
AltAlleleCount = data.frame(matrix(NA, nrow = nrow(vcf), ncol = 3))
colnames(AltAlleleCount) = c("SNP","Ref_counts","Alt_counts")
for (i in seq(1:nrow(vcf))){
    
  vcf_short[i,]=unfactor(vcf_short[i,])

    
  homozygote_alt_count <- length(which(vcf_short[i,]=="0|0"))
  homozygote_ref_count <- length(which(vcf_short[i,]=="1|1"))
  heterozygote_count <- length(which(vcf_short[i,]=="0|1" | vcf_short[i,]=="1|0"))

  #check
  check1 <- homozygote_alt_count + homozygote_ref_count+ heterozygote_count  == 685
  if (check1 != TRUE){
    print(paste("There is a problem at check1 at", i))

  }
    
    Alt_allelefreq <- round((homozygote_alt_count+heterozygote_count*0.5)/685,2)
    ref_allelefreq <- round((homozygote_ref_count+heterozygote_count*0.5)/685,2)

    check2 <- Alt_allelefreq+ref_allelefreq
    
    if (check2 !=1){
      print(paste("There is a problem at check2",i))
    }
    AltAlleleCount$SNP[i] = as.character(vcf[i,3])
    AltAlleleCount$Ref_counts[i] = homozygote_ref_count*2 + heterozygote_count
    AltAlleleCount$Alt_counts[i] = homozygote_alt_count*2 + heterozygote_count

  }
  
  
  write.table(AltAlleleCount,paste("AlternativeAlleleCount.txt"),col.names = T,row.names = F,quote=F,sep="\t")
  
  # plotting with non-polymorphic SNPs
  p1<-ggplot(AltAlleleCount, aes(x=Alt_counts)) + 
    geom_histogram(aes(y = ..count..),fill="black",binwidth=1,alpha=0.8)+
    xlab("Alternative allele count") +
    ylab("Number of SNPs") +
    theme_classic()
  p1
  filename_p1= paste("SFS_AllSNPs",Sys.Date(),".pdf",sep="")
  ggsave(filename_p1, plot = p1, width = 20, height = 14, unit = 'cm')
  
  # plotting only polymorphic SNPs
  rmv = which(AltAlleleCount$Alt_counts==0 | AltAlleleCount$Ref_counts==0)
  if (length(rmv)==0){
    AltAlleleCount_filt = AltAlleleCount
  } else{
    AltAlleleCount_filt = AltAlleleCount[-rmv,]
  }
  
  p15<-ggplot(AltAlleleCount_filt, aes(x=Alt_counts)) + 
    geom_histogram(aes(y = ..count..),fill="black",binwidth=1,alpha=0.8)+
    xlab("Alternative allele count") +
    ylab("Number of SNPs") +
    theme_classic()
  p15
  filename_p15= paste("SFS_AllPolymorphicSNPs_",Sys.Date(),".pdf",sep="")
  ggsave(filename_p15, plot = p15, width = 20, height = 14, unit = 'cm')
  

  # folded SFS
  AltAlleleCount_filt$MAC = AltAlleleCount_filt$Alt_counts
  for (i in seq(1:nrow(AltAlleleCount_filt))){
    AltAlleleCount_filt$MAC[i] = min(AltAlleleCount_filt[i,2:3])
  }
  p2<-ggplot(AltAlleleCount_filt, aes(x=MAC)) + 
    geom_histogram(aes(y = ..count..),fill="black",binwidth=1,alpha=0.8)+
    xlab("Alternative allele count") +
    ylab("Number of SNPs") +
    theme_classic()
  p2
  
  filename_p2= paste("FoldedSFS",Sys.Date(),".pdf",sep="")
  ggsave(filename_p2, plot = p2, width = 20, height = 14, unit = 'cm')






