# This should be a script to remove redundant samples within a panel
library(tidyverse)
library(stringr)
library(foreach)
args=commandArgs(trailingOnly = TRUE)



{
  genotypes=read.table(args[1], header = T, sep = "\t",comment.char = "?") #genotype file, not MAF filtered
  dim(genotypes)
  genotypes_shorter=genotypes[,10:ncol(genotypes)]
  genotypes_t=t(genotypes_shorter)
  dim(genotypes_t)
  rownames(genotypes_t)=colnames(genotypes_shorter)
  colnames(genotypes_t)=genotypes[,1]
  genotypes=genotypes_t
  names=rownames(genotypes)
}


# Look for redundant samples in general


# rewritten to take on row  (ind) at a time
Redundancy <- function(df,fixed_ind_number){
  m <- ncol(df) # number of loci
  n <- nrow(df)
  
  # fix one individual
  Fixed_Ind = df[fixed_ind_number,]
  currentIndName = rownames(df)[fixed_ind_number]
  print(paste("Now at individual:",currentIndName,sep=" "))
  
  results_matrix=data.frame(matrix(NA, nrow = n, ncol = 1))
  colnames(results_matrix)=currentIndName
  rownames(results_matrix)= rownames(df)
  
  for (j in seq(1:n)) {
    
    excluding <- which(Fixed_Ind=="./."| df[j,]=="./.")
    if (length(excluding)!=0){
      Eva = Fixed_Ind[-excluding]==df[j,-excluding]
      number_SNPs_being_comp = m-length(excluding)
    } else{
      Eva = Fixed_Ind==df[j,]
      number_SNPs_being_comp = m
    }
    number_complete_matches = sum(Eva)
    idx_to_check_for_incomplete_matches = which(Eva==FALSE)
    
    number_of_half_matches = 0
    
    if (length(idx_to_check_for_incomplete_matches)!=0){
      
      for (k in seq(1:length(idx_to_check_for_incomplete_matches))){
        genotype1 = Fixed_Ind[idx_to_check_for_incomplete_matches[k]]
        genotype2 = df[j,idx_to_check_for_incomplete_matches[k]]
        
        if (genotype1=="1|0" & genotype2=="0|1"){
          print(paste("A case where heterozygotes are written in two different ways,j=",j,"k=",k))
          number_of_half_matches = number_of_half_matches +2 #because it is a full match
        }
        
        if (genotype2=="1|0" & genotype1=="0|1"){
          print(paste("A case where heterozygotes are written in two different ways,j=",j,"k=",k))
          number_of_half_matches = number_of_half_matches +2 #because it is a full match
        }
        
        nuc1_geno1= str_split_fixed(genotype1, "|", 4)[2]
        nuc2_geno1= str_split_fixed(genotype1, "|", 4)[4]
        
        nuc1_geno2= str_split_fixed(genotype2, "|", 4)[2]
        nuc2_geno2= str_split_fixed(genotype2, "|", 4)[4]
        
        if (nuc1_geno1 ==nuc1_geno2){
          number_of_half_matches = number_of_half_matches+1
        }
        if (nuc2_geno1 ==nuc2_geno2){
          number_of_half_matches = number_of_half_matches+1
        }
      }
      
    }
    in_common <- number_complete_matches + 0.5*number_of_half_matches
    
    
    frac = in_common/number_SNPs_being_comp
    results_matrix[j,]=frac # return the fraction of the genotypes that two individuals have in common, half a match count 0.5 (that is "A/G" and "A/A" for instance), markers with missing entries are not included.
  }
  return(results_matrix)
}


library(doParallel)
myCluster <- makeCluster(8) 
registerDoParallel(myCluster)


Results = foreach(i=1:nrow(genotypes), .combine='cbind',.packages="stringr") %dopar% { 
  Redundancy(genotypes,i)
}
write.table(Results,file=paste("Redundancy_Matrix",Sys.Date(),sep="",".csv"),sep=",",col.names=T,row.names=T,quote=F)

