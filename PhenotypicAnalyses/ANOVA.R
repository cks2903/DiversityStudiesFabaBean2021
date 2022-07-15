# libraries ----
library(ggplot2)
library(tidyr)
library(lme4)
library(lmerTest)
library(harrypotter)
library(varhandle)



# read in data ----
file <- read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/_Phenotypes_qualityfiltered_minimum_2022-04-11_.csv",sep=",",head=T) #load the filtered phenotype file
dim(file)



#ANOVA FUNCTION ----
# ANOVA FUNCTION analyzing traits across different environments (defined as yearxloc combination AND yearxlocxdate combinations)
#be aware that every trait is treated as continous (this also accounts for hilum colour, seed colour etc. )
ANOVA <- function(file,trait){
  
  temp_df = file[which(file$PDID==trait),]
  number_of_environments = length(unique(temp_df$TRID))
  
  # dataframe for results ----
  results_df <- data.frame(trait = rep(NA,1),                     
                           Mean = rep(NA,1),
                           pvalue_mean = rep(NA,1),
                           VarG = rep(NA,1),
                           pvalue_G = rep(NA,1),
                           n_obs_G = rep(NA,1),
                           VarE = rep(NA,1),
                           pvalue_E = rep(NA,1),
                           n_obs_E = rep(NA,1),
                           VarGxE = rep(NA,1),
                           pvalue_GxE = rep(NA,1),
                           n_obs_GxE = rep(NA,1),
                           VarRep = rep(NA,1),
                           pvalue_Rep = rep(NA,1),
                           n_obs_Rep = rep(NA,1),
                           VarRes = rep(NA,1),
                           n_datapoints = rep(NA,1),
                           H2 = rep(NA,1)
  )
  
  # Remove NA scores
  
  if (length(which(is.na(temp_df$Score)==TRUE))!=0){
    temp_df <- temp_df[-which(is.na(temp_df$Score)==TRUE),]
  }
  
  # Make sure of factors
  temp_df$TRID <- as.factor(temp_df$TRID)
  temp_df$GPID <- as.factor(temp_df$GPID)
  temp_df$Rep <- as.factor(temp_df$Rep)
  
  
  if (number_of_environments!=1){
    # in this situation, we include effect of environment, GxE and nest Rep within and environment
      model <- lmer(Score ~ (1|GPID) + (1|TRID) + (1|GPID:TRID) + (1|Rep:TRID), data = temp_df) 
      estimates <- as.data.frame(VarCorr(model))
    
    # extract the relevant data and put into df
    VarGxE <- estimates$vcov[1]
    VarG <- estimates$vcov[2]
    VarRep <- estimates$vcov[3]
    VarE <- estimates$vcov[4]
    VarRes <- estimates$vcov[5]
  
    Mean <-  coef(summary(model))[1]
    pvalue_mean <-  coef(summary(model))[5]
    n_obs_GxE <- sapply(ranef(model),nrow)[1]
    n_obs_G <- sapply(ranef(model),nrow)[2]
    n_obs_Rep <- sapply(ranef(model),nrow)[3]
    n_obs_E <- sapply(ranef(model),nrow)[4]
    n_datapoints <- length(resid(model))
    
    

    
    
    # test effect of including single terms 
    model_refit <- lmer(Score ~ (1|GPID) + (1|TRID) + (1|GPID:TRID) + (1|Rep:TRID), data = temp_df,REML=F) 
  
    model_withoutGxE = lmer(Score ~ (1|GPID) + (1|TRID) + (1|Rep:TRID), data = temp_df,REML=F) 
    sig_GxE = anova(model_refit,model_withoutGxE)
    pvalue_GxE <- sig_GxE$`Pr(>Chisq)`[2]
  
    model_withoutRep = lmer(Score ~ (1|GPID) + (1|TRID) + (1|GPID:TRID), data = temp_df,REML=F) 
    sig_Rep = anova(model_refit,model_withoutRep)
    pvalue_Rep <- sig_Rep$`Pr(>Chisq)`[2]
  
    model_withoutG = lmer(Score ~  (1|Rep:TRID) + (1|TRID), data = temp_df,REML=F) 
    sig_G = anova(model_refit,model_withoutG)
    pvalue_G <- sig_G$`Pr(>Chisq)`[2]
  
    model_withoutE = lmer(Score ~ (1|GPID), data = temp_df,REML=F) 
    sig_E = anova(model_refit,model_withoutE)
    pvalue_E <- sig_E$`Pr(>Chisq)`[2]
    
    # calculate broad-sense heritability on a plot basis
    H2 = round(VarG/(VarG+VarGxE/n_obs_E+VarRes/n_obs_Rep),2)
}
  
  if (number_of_environments==1){
    # in this situation, we don't include any environment terms
      model <- lmer(Score ~ (1|GPID)  + (1|Rep), data = temp_df) 
    
    estimates <- as.data.frame(VarCorr(model))
    
    # extract the relevant data and put into df
    VarG <- estimates$vcov[1]
    VarGxE <- "NA"
    VarE <- "NA"
    pvalue_E <- "NA"
    pvalue_GxE <- "NA"
    VarRep <- estimates$vcov[2]
    VarRes <- estimates$vcov[3]
    
    Mean <-  coef(summary(model))[1]
    pvalue_mean <-  coef(summary(model))[5]
    n_obs_G <- sapply(ranef(model),nrow)[1]
    n_obs_GxE <- "NA"
    n_obs_E <- "NA"
    n_obs_Rep <- sapply(ranef(model),nrow)[2]
    n_datapoints <- length(resid(model))
    
      # test effect of including single terms
      model_refit = lmer(Score ~ (1|GPID) + (1|Rep), data = temp_df,REML=F)
      model_withoutRep = lmer(Score ~ (1|GPID), data = temp_df,REML=F) 
      sig_Rep = anova(model_refit,model_withoutRep)
      pvalue_Rep <- sig_Rep$`Pr(>Chisq)`[2]
      
      model_withoutG = lmer(Score ~  (1|Rep), data = temp_df,REML=F) 
      sig_G = anova(model_refit,model_withoutG)
      pvalue_G <- sig_G$`Pr(>Chisq)`[2]
      
      sig_E = "NA"
      sig_GxE ="NA"
    
    # calculate broad-sense heritability on a plot basis
    H2 = round(VarG/(VarG+VarRes/n_obs_Rep),2)
}
  
  
  # output should be a row for a dataframe with all necessary information
  
  results_df$trait = trait
  results_df$Mean = Mean
  results_df$pvalue_mean = pvalue_mean
  results_df$VarG = VarG
  results_df$pvalue_G = pvalue_G
  results_df$n_obs_G = n_obs_G
  results_df$VarE = VarE
  results_df$pvalue_E = pvalue_E
  results_df$n_obs_E = n_obs_E
  results_df$VarGxE = VarGxE
  results_df$pvalue_GxE = pvalue_GxE
  results_df$n_obs_GxE = n_obs_GxE
  results_df$VarRep = VarRep
  results_df$pvalue_Rep = pvalue_Rep
  results_df$n_obs_Rep = n_obs_Rep
  results_df$VarRes = VarRes
  results_df$n_datapoints = n_datapoints
  results_df$H2 =H2
  
  return(results_df)
}


# Run multi-environmental analyses for all traits ----
list_of_traits_filtered <- unique(file$PDID)

for (i in seq(1:length(list_of_traits_filtered))){
  
  # if the merged dataset doesn't exist, create it
    if (!exists("All_result")){
      All_result <- ANOVA(file,list_of_traits_filtered[i])
    }
    
    # if the merged dataset does exist, append to it
    if (exists("All_result")){
      temp_results <-ANOVA(file,list_of_traits_filtered[i])
      All_result<-rbind(All_result, temp_results)
      rm(temp_results)
    }
}



# Create plots showing variance for each trait ----
All_result
All_result_copy = All_result
All_result_copy$VarE = as.numeric(as.character(All_result_copy$VarE))
All_result_copy$VarG = as.numeric(as.character(All_result_copy$VarG))
All_result_copy$VarGxE = as.numeric(as.character(All_result_copy$VarGxE))
All_result_copy$VarRep= as.numeric(as.character(All_result_copy$VarRep))
All_result_copy$VarRes= as.numeric(as.character(All_result_copy$VarRes))

All_result_Ordered_after_H2 <- All_result_copy[order(-All_result_copy$H2),]
All_result_Ordered_after_H2=All_result_Ordered_after_H2[-19,]

LongerDF=pivot_longer(All_result_Ordered_after_H2, c("VarG","VarE","VarGxE","VarRep","VarRes"),names_to="variance",values_to="Score")


# Stacked + percent
# change order of levels 
LongerDF$variance <- factor(LongerDF$variance,levels=c("VarRes","VarRep","VarGxE","VarE","VarG"))
LongerDF$trait <- factor(LongerDF$trait , levels = All_result_Ordered_after_H2$trait)



All_result_Ordered_after_H2 <- LongerDF[order(-LongerDF$H2),]
guide <- All_result_Ordered_after_H2$trait[!duplicated(All_result_Ordered_after_H2$trait)]

LongerDF$trait <- factor(LongerDF$trait , levels = guide)
LongerDF_oneEnviron = LongerDF[which(LongerDF$n_obs_E=="NA"),]
LongerDF_MoreEnviron = LongerDF[-which(LongerDF$n_obs_E=="NA"),]

# split into two dataframes depending on whether or not multiple environments were included
#for paper

plot1 <- ggplot(LongerDF_MoreEnviron, aes(fill=variance, y=Score, x=trait)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_hp(discrete = TRUE, option = "Ravenclaw", name = "Cut")  +
  theme_classic() +  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  xlab("Trait") +
  ylab("Fraction of variance explained")

filename <-paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/ANOVA/MultiEnvironmentalTraits",Sys.Date(),".pdf",sep="_")
ggsave(filename,plot1,width=20,height=12,units="cm")


LongerDF_oneEnviron$variance = unfactor(LongerDF_oneEnviron$variance)
LongerDF_oneEnviron = LongerDF_oneEnviron[-which(LongerDF_oneEnviron$variance=="VarGxE" | LongerDF_oneEnviron$variance=="VarE" ),]
LongerDF_oneEnviron$variance =factor(LongerDF_oneEnviron$variance,levels=c("VarRes","VarRep","VarG"))


plot1.5 <- ggplot(LongerDF_oneEnviron, aes(fill=variance, y=Score, x=trait)) + 
  geom_bar(position="fill", stat="identity")+ scale_fill_manual(values = c("VarRes" = "#005B8E",
                                                                           "VarRep" = "#2988BE",
                                                                           "VarG" = "#AA4F00"))  +
  theme_classic() +  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  xlab("Trait") +
  ylab("Fraction of variance explained")

filename <-paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/ANOVA/OneTrial",Sys.Date(),".pdf",sep="_")
ggsave(filename,plot1.5,width=20,height=12,units="cm")

# 

# Save the analysis ----
filename2 <-paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/ANOVA/","Variances_MULTIENVIRONMENTAL",Sys.Date(),".txt",sep="_")
write.table(All_result,filename2,col.names = T,row.names = F,quote=F,sep="\t")



# Run one trait, one trial each time ----
file$TRID_PDID = paste(file$TRID,file$PDID,sep=":")
unique(file$TRID_PDID ) # 71 unique trait and trial combos
filtered_file10 = file
filtered_file10$TRID = filtered_file10$TRID_PDID
filtered_file10$PDID = filtered_file10$TRID_PDID

All_result =NULL
list_of_traits <- unique(filtered_file10$TRID_PDID)

for (i in seq(1:length(list_of_traits))){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("All_result")){
    All_result <- ANOVA(filtered_file10,list_of_traits[i])
  }
  
  # if the merged dataset does exist, append to it
  if (exists("All_result")){
    temp_results <-ANOVA(filtered_file10,list_of_traits[i])
    All_result<-rbind(All_result, temp_results)
    rm(temp_results)
  }
}



# Create plots showing variance for each trait ----
All_result_copy = All_result
All_result_copy$VarE = as.numeric(as.character(All_result_copy$VarE))
All_result_copy$VarG = as.numeric(as.character(All_result_copy$VarG))
All_result_copy$VarGxE = as.numeric(as.character(All_result_copy$VarGxE))
All_result_copy$VarRep= as.numeric(as.character(All_result_copy$VarRep))
All_result_copy$VarRes= as.numeric(as.character(All_result_copy$VarRes))

All_result_Ordered_after_H2 <- All_result_copy[order(-All_result_copy$H2),]


LongerDF=pivot_longer(All_result_Ordered_after_H2, c("VarG","VarE","VarGxE","VarRep","VarRes"),names_to="variance",values_to="Score")


# Stacked + percent
# change order of levels 
LongerDF$variance <- factor(LongerDF$variance,levels=c("VarRes","VarRep","VarGxE","VarE","VarG"))
LongerDF$trait <- factor(LongerDF$trait , levels = All_result_Ordered_after_H2$trait)

# make better names for traits
LongerDF$trait <- unfactor(LongerDF$trait)


for (i in (1:nrow(LongerDF))){
  trial <- strsplit(LongerDF$trait[i], ":")[[1]][1]
  trait <- strsplit(LongerDF$trait[i], ":")[[1]][2]
  
  if (trial=="29"){
    normal_trial_name = "Sejet 2021"
  }
  if (trial=="28"){
    normal_trial_name = "Dyngby 2020"
  }
  
  if (trial=="27"){
    normal_trial_name = "Sejet 2020"
  }
  
  if (trial=="31"){
    normal_trial_name = "Sejet 2020"
  }
  
  normal_trait_name = dic$trait[which(dic$PDID==trait)]
  LongerDF$trait[i] = paste(normal_trait_name, normal_trial_name,sep=" ")
  
}


All_result_Ordered_after_H2 <- LongerDF[order(-LongerDF$H2),]
guide <- All_result_Ordered_after_H2$trait[!duplicated(All_result_Ordered_after_H2$trait)]

LongerDF$trait <- factor(LongerDF$trait , levels = guide)


# Save the analysis ----
filename20 <-paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/ANOVA/","Variances_singlenvironments_filteredforbadtrials",Sys.Date(),".txt",sep="_")
write.table(All_result,filename20,col.names = T,row.names = F,quote=F,sep="\t")










