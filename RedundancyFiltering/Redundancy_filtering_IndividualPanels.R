# This is a script to remove redundant samples within a panel
library(tidyverse)
library(stringr)
library(foreach)
args=commandArgs(trailingOnly = TRUE)


{
  genotypes=read.table(args[1], header = T, sep = "\t",comment.char = "#") #genotype file, not MAF filtered
  dim(genotypes)
  genotypes_shorter=genotypes[,10:ncol(genotypes)]
  genotypes_t=t(genotypes_shorter)
  dim(genotypes_t)
  rownames(genotypes_t)=colnames(genotypes_shorter)
  colnames(genotypes_t)=genotypes[,1]
  genotypes=genotypes_t
  names=rownames(genotypes)
}


# divide df for each panel
Dic <- read.table("Lines_Panel_DicNew.txt",head=T,sep="\t") # A file that indicates which panel a given sample belong to


Panels = unique(Dic[,2])

Accessions_EUCLEG = Dic[which(Dic[,2]=="EUCLEG"),1]
Accessions_RSBP = Dic[which(Dic[,2]=="RSBP"),1]
Accessions_KHAZAEI_4PRIL = Dic[which(Dic[,2]=="Four-way-cross"),1]
Accessions_NorfabCore = Dic[which(Dic[,2]=="NORFAB"),1]
Accessions_GWB = Dic[which(Dic[,2]=="GWB"),1]
Accessions_PROFABA = Dic[which(Dic[,2]=="ProFaba"),1]
Accessions_MAGIC = Dic[which(Dic[,2]=="Seven-paren-MAGIC"),1]
Accessions_VICCI = Dic[which(Dic[,2]=="VICCI"),1]


df_EUCLEG <- genotypes[which(rownames(genotypes) %in% Accessions_EUCLEG),]
df_RSBP <- genotypes[which(rownames(genotypes) %in% Accessions_RSBP),]
df_KHAZAEI_4PRIL <- genotypes[which(rownames(genotypes) %in% Accessions_KHAZAEI_4PRIL),]
df_NorfabCore <- genotypes[which(rownames(genotypes) %in% Accessions_NorfabCore),]
df_GWB <- genotypes[which(rownames(genotypes) %in% Accessions_GWB),]
df_PROFABA<- genotypes[which(rownames(genotypes) %in% Accessions_PROFABA),]
df_MAGIC<- genotypes[which(rownames(genotypes) %in% Accessions_MAGIC),]
df_VICCI<- genotypes[which(rownames(genotypes) %in% Accessions_VICCI),]


Redundancy_within_panel <- function(panel_df,fixed_ind_number){
  m <- ncol(panel_df) # number of loci
  n <- nrow(panel_df)
  
  # fix one individual
  Fixed_Ind = panel_df[fixed_ind_number,]
  currentIndName = rownames(panel_df)[fixed_ind_number]
  print(paste("Now at individual:",currentIndName,sep=" "))
  
  results_matrix=data.frame(matrix(NA, nrow = n, ncol = 1))
  colnames(results_matrix)=currentIndName
  rownames(results_matrix)= rownames(panel_df)
    
    for (j in seq(1:n)) {
      
      excluding <- which(Fixed_Ind=="./."| panel_df[j,]=="./.")
      if (length(excluding)!=0){
        Eva = Fixed_Ind[-excluding]==panel_df[j,-excluding]
        number_SNPs_being_comp = m-length(excluding)
      } else{
        Eva = Fixed_Ind==panel_df[j,]
        number_SNPs_being_comp = m
      }
      number_complete_matches = sum(Eva)
      idx_to_check_for_incomplete_matches = which(Eva==FALSE)
      
      number_of_half_matches = 0
      
      if (length(idx_to_check_for_incomplete_matches)!=0){
        
        for (k in seq(1:length(idx_to_check_for_incomplete_matches))){
          genotype1 = Fixed_Ind[idx_to_check_for_incomplete_matches[k]]
          genotype2 = panel_df[j,idx_to_check_for_incomplete_matches[k]]
          
          nuc1_geno1= str_split_fixed(genotype1, "/", 2)[,1]
          nuc2_geno1= str_split_fixed(genotype1, "/", 2)[,2]
          
          nuc1_geno2= str_split_fixed(genotype2, "/", 2)[,1]
          nuc2_geno2= str_split_fixed(genotype2, "/", 2)[,2]
          
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



# run in parallel style
library(doParallel)
myCluster <- makeCluster(8) 
registerDoParallel(myCluster)


# tell the script which dataframe is wanted (args[2]) for instance df_EUCLEG
wanted_df = get(args[2])
Results = foreach(i=1:nrow(wanted_df), .combine='cbind',.packages="stringr") %dopar% { 
  Redundancy_within_panel(wanted_df,i)
}
write.table(Results,file=paste("Redundancy_Matrix",args[1],Sys.Date(),sep="",".csv"),sep=",",col.names=T,row.names=T,quote=F)

