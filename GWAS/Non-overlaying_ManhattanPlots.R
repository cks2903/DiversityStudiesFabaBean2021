#########################################
#########################################
# Make a Manhattan plot of any .pvals file
# inputs are:
# 1) .pvals file
# 2) name of trait
#########################################
#########################################

library(ggplot2)
library(dplyr)
args=commandArgs(trailingOnly = TRUE)

# load file

data <- read.table(args[1], header=TRUE, sep=",")

datafiltered1=data

trait=args[2]

# SimpleM threshold ---
library(matrixStats)

#============================================================================
# simpleM_Ex.R 

#============================================================================
# License:  GPL version 2 or newer. 
# NO warranty. 

#============================================================================
# citation: 
#
# Gao X, Starmer J and Martin ER (2008) A Multiple Testing Correction Method for
# Genetic Association Studies Using Correlated Single Nucleotide Polymorphisms. 
# Genetic Epidemiology 32:361-369
#
# Gao X, Becker LC, Becker DM, Starmer J, Province MA (2009) Avoiding the high 
# Bonferroni penalty in genome-wide association studies. Genetic Epidemiology 
# (Epub ahead of print) 

#============================================================================
# readme: 
# example SNP file format:
# row => SNPs
# column => Unrelated individuals 

# The data file should contain only POLYMORPHIC SNPs. 

# Missing values should be imputed. 
# There should be NO missing values in the SNP data file.
# SNPs are coded as 0, 1 and 2 for the number of reference alleles. 
# SNPs are separated by one-character spaces. 

# You may need to change file path (search for "fn_In" variable) 
# depending on where your snp file is stored at.

#============================================================================
# Meff through the PCA approach 
# use a part of the eigen values according to how much percent they contribute
# to the total variation 
Meff_PCA <- function(eigenValues, percentCut){
  totalEigenValues <- sum(eigenValues)
  myCut <- percentCut*totalEigenValues
  num_Eigens <- length(eigenValues)
  myEigenSum <- 0
  index_Eigen <- 0
  
  for(i in 1:num_Eigens){
    if(myEigenSum <= myCut){
      myEigenSum <- myEigenSum + eigenValues[i]
      index_Eigen <- i
    }
    else{
      break
    }
  }	
  return(index_Eigen)
}

#============================================================================
# infer the cutoff => Meff
inferCutoff <- function(dt_My){
  CLD <- cor(dt_My)
  eigen_My <- eigen(CLD)
  
  # PCA approach
  eigenValues_dt <- abs(eigen_My$values)
  Meff_PCA_gao <- Meff_PCA(eigenValues_dt, PCA_cutoff)
  return(Meff_PCA_gao)
}

#============================================================================
PCA_cutoff <- 0.995

#============================================================================
# fix length, simpleM
fn_In <- "../Genotypes_FarmCPU.txt"				# <---- change path here!!!
cls <- c("character",rep("numeric",21345))
mySNP_nonmissing <- read.table(fn_In,sep="\t",stringsAsFactors = F,colClasses=cls,head=T,check.names=F)	
mySNP_nonmissing <- mySNP_nonmissing[,2:ncol(mySNP_nonmissing)]




# remove SNPs not used because MAF < 0.05
SNPinfo <- read.table("../GenotypeInfo_FarmCPU.txt",sep="\t",stringsAsFactors = F,head=T) 
KEEP = SNPinfo[which(SNPinfo[,1] %in% datafiltered1$SNP),1]
Keepidx = which(colnames(mySNP_nonmissing) %in% KEEP)
mySNP_nonmissing = mySNP_nonmissing[,Keepidx]
mySNP_nonmissing <- t(mySNP_nonmissing)



# remove monomorphic SNPs
mySNP_nonmissing= mySNP_nonmissing[-which(rowVars(as.matrix(mySNP_nonmissing))==0),]


numLoci <- length(mySNP_nonmissing[, 1])

simpleMeff <- NULL
fixLength <- 133 
i <- 1
myStart <- 1
myStop <- 1
while(myStop < numLoci){
  myDiff <- numLoci - myStop 
  if(myDiff <= fixLength) break
  
  myStop <- myStart + i*fixLength - 1
  snpInBlk <- t(mySNP_nonmissing[myStart:myStop, ])
  MeffBlk <- inferCutoff(snpInBlk)
  simpleMeff <- c(simpleMeff, MeffBlk)
  myStart <- myStop+1
}
snpInBlk <- t(mySNP_nonmissing[myStart:numLoci, ])
MeffBlk <- inferCutoff(snpInBlk)
simpleMeff <- c(simpleMeff, MeffBlk)

cat("Total number of SNPs is: ", numLoci, "\n")
cat("Inferred Meff is: ", sum(simpleMeff), "\n")

effective_number_of_loci <- sum(simpleMeff)
#============================================================================
# end 




# Plotting
{
  #prepare for plotting
  datafiltered <- subset(datafiltered1, Chromosome >0)

  datafiltered$p.bh <- p.adjust(datafiltered$P.value, "BH")
  p.bonferroni <- 1/nrow(datafiltered)*0.05
  p.effectivenumberofloci <- 1/effective_number_of_loci*0.05
  datafiltered=datafiltered[which(datafiltered$maf>=0.01),]
 
  
  
  don <- datafiltered %>% 
    
    # Compute chromosome size
    group_by(Chromosome) %>% 
    summarise(chr_len=max(Position)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(datafiltered, ., by=c("Chromosome"="Chromosome")) %>%
    
    # Add a cumulative position of each SNP
    arrange(Chromosome, Position) %>%
    mutate( BPcum=Position+tot) %>%
    
    mutate( is_highlight=ifelse(p.bh<0.05, "yes", "no")) 
  
  
  # Prepare X axis 
  axisdf = don %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  
  # now calculate genomic inflation factor
chisq <- qchisq(1-don$P.value,1)
lampda=median(chisq)/qchisq(0.5,1) # Genomic inflation factor is 0.99 <1. No signs of systematic bias
print(paste("Genomic inflation factor of", trait, "is:", round(lampda,2),sep=" "))

if (lampda<0.9 | lampda>1.1){
  don$P.value = don$P.value/lampda
}

  
  # plotting
  p <- ggplot(don, aes(x=BPcum, y=-log10(P.value))) +
    
    # Show all points
    geom_point( aes(color=as.factor(Chromosome)), alpha=0.8, size=2) +
    scale_color_manual(values = rep(c("gray33", "gray67"), 16 )) +
    
    # custom X axis:
    scale_x_continuous(label = c("1S","1L","2","3","4","5","6"), breaks= axisdf$center, expand = c(0, 0)) +
    #scale_y_continuous(limits = c(0,7.5), expand = expand_scale(mult = c(0, .1))) +     # remove space between plot area and x axis
    # Add thresholds
    geom_hline(yintercept=-log10(p.effectivenumberofloci), linetype="dashed", color = "green",size=1.5) + 
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
     panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 11),
      title = element_text(size = 10)
    ) +
    labs(y = ("-log10(p-value)"), x = ("Chromosome"), title = (trait))
  
  p
}
ggsave(paste('ManhattanPlot_',trait,Sys.Date(),'.pdf',sep="_"), plot = p, width = 30, height = 7.5, unit = 'cm')
ggsave(paste('ManhattanPlot_',trait,Sys.Date(),'.png',sep="_"), plot = p, width = 30, height = 7.5, unit = 'cm')





