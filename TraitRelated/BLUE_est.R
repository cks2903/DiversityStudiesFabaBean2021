# libraries ----
#library(ggplot2)
#library(tidyr)
library(lme4)
library(nlme)
#library(lmerTest)
#library(harrypotter)
#library(varhandle)


# read in data ----

# the full phenotype file that has already been filtered for unusable data
file <- read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/ANOVA/_FilteredPhenotypefile_2022-04-11_.txt",sep=";",head=T) #phenotype scores

# the file with information about whether effect are significant in a multienvironmental field trial, generated in the ANOVA analysis
sigfile_multi <- read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/ANOVA/_Variances_MULTIENVIRONMENTAL_filteredforbadtrials_2022-03-25_.txt",sep="\t",head=T)



for (i in (1:nrow(file))){
  trial <- strsplit(file$TRID_PDID[i], ":")[[1]][1]
  trait <- strsplit(file$TRID_PDID[i], ":")[[1]][2]
  
  if (trial=="29"){
    normal_trial_name = "Se_2021"
  }
  if (trial=="28"){
    normal_trial_name = "NS_2020"
  }
  
  if (trial=="28.1"){
    normal_trial_name = "NS_2020_early_date"
  }
  
  if (trial=="28.2"){
    normal_trial_name = "NS_2020_late_date"
  }
  
  if (trial=="27"){
    normal_trial_name = "Se_2020"
  }
  
  if (trial=="31"){
    normal_trial_name = "NS_2021"
  }
  
  normal_trait_name = dic$trait[which(dic$PDID==trait)]
  file$TRID_PDID[i] = paste(normal_trait_name, normal_trial_name,sep="_")
  
}





# Phenotype averages pr. trial:trait combo ---- 
Phen_avg <- function(File,trait){
  
  current_df = File[which(File$TRID_PDID==trait),]
  n = nrow(current_df)
  
  Resulting_table <- aggregate(current_df$Score, list(current_df$Name), mean,na.action=na.omit)
  sec_col_name<- paste(trait,"phen_avg",sep="_")
  colnames(Resulting_table) = c("Line",sec_col_name)
  
  return(Resulting_table)
}

# apply above function ----
Phen_avg_df <- file
list_of_traits_1 <- (unique(Phen_avg_df$TRID_PDID)) 

for (i in seq(1:length(list_of_traits_1))){
    
# if the merged dataset doesn't exist, create it
  if (i==1){
    All_result <- Phen_avg(Phen_avg_df,list_of_traits_1[i])
  }
    
  # if the merged dataset does exist, append to it
  if (i!=1){
    temp_results <-Phen_avg(Phen_avg_df,list_of_traits_1[i])
    All_result<-merge(All_result, temp_results,by="Line",all=T)
    
    rm(temp_results)
  }
}
  
All_result_sorted <- All_result[order(All_result$Line),]

filename1 <- paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/BLUEs_and_more_phenotypes/","AveragesPrEnvironment",Sys.Date(),".csv",sep="_")
write.table(All_result_sorted,filename1,col.names = TRUE, row.names = F, quote=F,sep=",")


# BLUEs multienvironmental trials function ----
BLUEs_multienvironmental <- function(File,sigfile_multi,trait){
  
    TEMP_DF_obs = File[which(File$PDID==trait),]
    
    Relevance_of_Env <- sigfile_multi$pvalue_E[which(sigfile_multi$trait==trait)]
    Relevance_of_GxE <- sigfile_multi$pvalue_GxE[which(sigfile_multi$trait==trait)]
    Relevance_of_Rep <- sigfile_multi$pvalue_Rep[which(sigfile_multi$trait==trait)]
    
    TEMP_DF_obs$Rep = as.factor( TEMP_DF_obs$Re)
    TEMP_DF_obs$Name = as.factor(TEMP_DF_obs$Name)
    TEMP_DF_obs$TRID = as.factor(TEMP_DF_obs$TRID)
    
    # model for relevance of environments, GxE and reps
    if (Relevance_of_Env< 0.05 & Relevance_of_GxE< 0.05 &  Relevance_of_Rep<0.05){
          
      model <- lmer(Score ~ Name + (1|TRID) + (1|Name:TRID) + (1|Rep:TRID), data = TEMP_DF_obs) 
      BLUEs <- fixef(model)[2:length( fixef(model))]
    }
      
    # model for relevance of environments, GxE
    if (Relevance_of_Env< 0.05 & Relevance_of_GxE< 0.05 &  Relevance_of_Rep>=0.05){
        
      model <- lmer(Score ~ Name + (1|TRID) + (1|Name:TRID), data = TEMP_DF_obs) 
      BLUEs <- fixef(model)[2:length( fixef(model))]
    }
      
    # model for relevance of environments, reps
    if (Relevance_of_Env< 0.05 & Relevance_of_GxE>=0.05 &  Relevance_of_Rep<0.05){
        
      model <- lmer(Score ~ Name + (1|TRID) + (1|Rep:TRID), data = TEMP_DF_obs) 
      BLUEs <- fixef(model)[2:length( fixef(model))]
    }
      
    # model for relevance of GxE and reps
    if (Relevance_of_Env>= 0.05 & Relevance_of_GxE< 0.05 &  Relevance_of_Rep<0.05){
        
      model <- lmer(Score ~ Name + (1|Name:TRID) + (1|Rep:TRID), data = TEMP_DF_obs) 
      BLUEs <- fixef(model)[2:length( fixef(model))]
      }
      
    # model for relevance of E
    if (Relevance_of_Env< 0.05 & Relevance_of_GxE>=0.05 &  Relevance_of_Rep>=0.05){
        
      model <- lmer(Score ~ Name + (1|TRID), data = TEMP_DF_obs) 
      BLUEs <- fixef(model)[2:length( fixef(model))]
    }
      
    # model for relevance of GxE
    if (Relevance_of_Env>= 0.05 & Relevance_of_GxE<0.05 &  Relevance_of_Rep>=0.05){
        
      model <- lmer(Score ~ Name + (1|Name:TRID), data = TEMP_DF_obs) 
      BLUEs <- fixef(model)[2:length( fixef(model))]
    }
      
    # model for relevance of Reps
    if (Relevance_of_Env>= 0.05 & Relevance_of_GxE>=0.05 &  Relevance_of_Rep<0.05){
      model <- lmer(Score ~ Name + (1|Rep), data = TEMP_DF_obs) 
      BLUEs <- fixef(model)[2:length( fixef(model))]
    }
      
    # model for relevance of nothing
    if (Relevance_of_Env>= 0.05 & Relevance_of_GxE>=0.05 &  Relevance_of_Rep>=0.05){
      model <- gls(Score ~ Name, data = TEMP_DF_obs) 
      BLUEs <- model$coefficients[2:length(model$coefficients)]
    }
      
    BLUE_df = as.data.frame(BLUEs)
    BLUE_df2 = cbind(rownames(BLUE_df),BLUE_df$BLUEs)
    col2name = paste("BLUE_trait",trait,sep="_")
    colnames(BLUE_df2)= c("Line",col2name)
    return(BLUE_df2)
      
    }


# Apply function above ----
# Multi-environmental trials_BLUEs
BLUEs_MultiFile <- file
listoftraits_multienvironmental <- unique(file$PDID)

# remove the traits only present in one environment
rmv1 = which(is.na(sigfile_multi$VarE)=="TRUE") # only present in one environment

listoftraits_multienvironmental_filt = listoftraits_multienvironmental[-c(rmv1)] 


for (i in seq(1:length(listoftraits_multienvironmental_filt))){
  
  # if the merged dataset doesn't exist, create it
  if (i==1){
    All_result_multi <- BLUEs_multienvironmental(BLUEs_MultiFile,sigfile_multi,listoftraits_multienvironmental_filt[i])
  }
  
  # if the merged dataset does exist, append to it
  if (i!=1){
    temp_results <-BLUEs_multienvironmental(BLUEs_MultiFile,sigfile_multi,listoftraits_multienvironmental_filt[i])
    All_result_multi<-merge(All_result_multi, temp_results,by="Line",all=T)
    
    rm(temp_results)
  }
}


All_result_multi$Line <- substring(All_result_multi$Line, 5)      # remove name part of line name
All_result_multi_sorted <- All_result_multi[order(All_result_multi$Line),]

filename2 <- paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/BLUEs_and_more_phenotypes/","BLUEs_multienvironment",Sys.Date(),".csv",sep="_")
write.table(All_result_multi_sorted,filename2,col.names = TRUE, row.names = F, quote=F,sep=",")






# BLUEs for each trail:trial combo ----

BLUEs_eachTRIDPDIDcombo <- file

# This only makes sense to do if REP has an influence on the trait
ANOVA_singletrials = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/ANOVA/_Variances_singlenvironments_filteredforbadtrials_2022-03-25_.txt",sep="\t",head=T)

listoftraits_TRID_PDIDcombos_filtered = ANOVA_singletrials$trait[which(ANOVA_singletrials$pvalue_Rep<0.05)]


File=BLUEs_eachTRIDPDIDcombo
sigfile_multi = ANOVA_singletrials
trait = listoftraits_TRID_PDIDcombos_filtered[1]


# BLUEs multienvironmental trialsfunction ----
# Function to refit and extract BLUEs from traits scored in more environments 
BLUEs_singleenv<- function(File,sigfile_multi,trait){
  
  TEMP_DF_obs = File[which(File$TRID_PDID==trait),]
  
  Relevance_of_Rep <- sigfile_multi$pvalue_Rep[which(sigfile_multi$trait==trait)]
  
  TEMP_DF_obs$Rep = as.factor( TEMP_DF_obs$Rep)
  TEMP_DF_obs$Name = as.factor(TEMP_DF_obs$Name)
  
  # model for relevance of environments, GxE and reps
  if (Relevance_of_Rep< 0.05){
    
    model <- lmer(Score ~ Name + (1|Rep), data = TEMP_DF_obs) 
    BLUEs <- fixef(model)[2:length( fixef(model))]
  } else {
    print("this trait do not  have significant rep variation")
  }
  
  BLUE_df = as.data.frame(BLUEs)
  BLUE_df2 = cbind(rownames(BLUE_df),BLUE_df$BLUEs)
  col2name = paste("BLUE_trait",trait,sep="_")
  colnames(BLUE_df2)= c("Line",col2name)
  return(BLUE_df2)
  
}

# run function

for (i in seq(1:length(listoftraits_TRID_PDIDcombos_filtered))){
  
  # if the merged dataset doesn't exist, create it
  if (i==1){
    All_result_singleenv <- BLUEs_singleenv(BLUEs_eachTRIDPDIDcombo,ANOVA_singletrials,listoftraits_TRID_PDIDcombos_filtered[i])
  }
  
  # if the merged dataset does exist, append to it
  if (i!=1){
    temp_results <-BLUEs_singleenv(BLUEs_eachTRIDPDIDcombo,ANOVA_singletrials,listoftraits_TRID_PDIDcombos_filtered[i])
    All_result_singleenv<-merge(All_result_singleenv, temp_results,by="Line",all=T)
    
    rm(temp_results)
  }
}

head(All_result_singleenv)
dim(All_result_singleenv)
All_result_singleenv$Line <- substring(All_result_singleenv$Line, 5)      # remove name part of line name
head(All_result_singleenv)

All_result_singleenv_sorted <- All_result_singleenv[order(All_result_singleenv$Line),]
head(All_result_singleenv_sorted)

filename3 <- paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/BLUEs_and_more_phenotypes/","BLUEs_singleenvironment",Sys.Date(),".csv",sep="_")
write.table(All_result_singleenv_sorted,filename3,col.names = TRUE, row.names = F, quote=F,sep=",")




