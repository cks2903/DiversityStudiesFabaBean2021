
#libraries ----
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(ggpubr)
library("patchwork")



# Load data ----
df <- read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Se20_NS20_Se21_NS21.csv",sep=",",head=T) # raw genotype data from field trials
dim(df)
head(df)

dic <- read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/binwidth_dic.txt",head=T,sep="\t")

# Initial filtering ----

# start by removing NA scores
rmv1 <- which(is.na(df$Score)==TRUE | df$Score=="")
df_filt1 <- df[-rmv1,]

# Count number of unique lines
uniquelines <- unique(df_filt1$Name)
Number_of_uniquelines <- length(uniquelines) 

MAGIClines <- uniquelines[grepl( "MAGIC", uniquelines, fixed = TRUE)]
MAGIClines_n <- length(MAGIClines) 

CommercialLines <- uniquelines[which(grepl( "MAGIC", uniquelines, fixed = TRUE)==FALSE)]
Commerciallines_n <- length(CommercialLines) # 10-11 repeated commercial check cultivars

df_filt1 = df_filt1[which(df_filt1$Name %in% MAGIClines),] #remove the commercial lines
# How many of the lines are in all 3 trials

count_all_trials <-  0
count_three_trials <-  0
count_two_trials <- 0
count_one_trials <- 0
for (i in seq(1:Number_of_uniquelines)){
  temp_df <- df_filt1[which(df_filt1$Name == uniquelines[i]),]
  trials_with_line_ = unique(temp_df$TRID)
  trials_with_line <- length(unique(round(trials_with_line_))) # so that same line scored two different date in a trial is not counted twice
  
  if (trials_with_line ==4){
    count_all_trials = count_all_trials+1
  }
  
  if (trials_with_line ==3){
    count_three_trials = count_all_trials+1
  }
  
  if (trials_with_line ==2){
    count_two_trials = count_two_trials+1
  }
  
  if (trials_with_line ==1){
    count_one_trials = count_one_trials+1
  }
  
}
count_all_trials #207 lines are in all 4 trials
count_three_trials #0 lines in precisely 3 trials
count_two_trials #41 lines are in precisely 2 trials
count_one_trials #2 lines are in precisely 1 trial only


# Initial filtering ----
# Remove accessions where we have less than two observations for a trait
df_filt1$PDID_Name =paste(df_filt1$PDID,df_filt1$Name,sep=":")

rmv2 <- as.data.frame(which(table(df_filt1$PDID_Name)<2))
rmv2_ = row.names(rmv2)

if (length(rmv2_)!=0){
  rmv2__ <- which(df_filt1$PDID_Name %in% rmv2_)
  df_filt2 <- df_filt1[-rmv2__,] # removed 6 observations
}


# function that makes histograms
# Does this by trait and by trial
histogrammer <-function(data,path){
  
  # filter dataframe to include only one trait
  
  traits_Scored <- unique(data$PDID)
  traits_n <- length(traits_Scored)
  

  for (trait in traits_Scored){

    trait_df = data[which(data$PDID==trait),]
    trait_df$Score = as.numeric(as.character(trait_df$Score))
    
    if (length(which(is.na(trait_df$Score)==T))!=0){
      trait_df = trait_df[-which(is.na(trait_df$Score)==T),]
    }
    Unique_trials <- unique(trait_df$TRID)
    Number_of_trials <- length(Unique_trials)
    
    # colours for different trial
    colour_df =as.data.frame(cbind(c(28.1,28.2,28.0,27.0,29.,31.0),c("#10108B","#10108B","#10108B","#104E8B","#108B8B","pink")))
    
    
    # this loops create a unique df for each TRID
    for (j in seq(1:Number_of_trials)){
      Current_TRID <- Unique_trials[j]
      df_x <- trait_df[which(trait_df$TRID==Current_TRID),]
      assign(paste("df",j,sep=""),df_x)
    }
    
    # plotting
    
    #min = dic[which(dic$PDID==trait),3]
    #max = dic[which(dic$PDID==trait),4]
    min=NULL
    max=NULL
    
    if (length(min)==0){
      min = min(trait_df$Score)
    }
    
    if (length(max)==0){
      max = max(trait_df$Score)
    }
    
    
    if (length(which(trait_df$Score==99))!=0){
      max = 99
    }
    
    if (is.null(min) | is.na(min)){
      min = min(trait_df$Score)
    }
    
    if (is.null(max)| is.na(max)){
      max = max(trait_df$Score)
    }
 
    
    # define binwidth based on PDID
    
    binwidth_custom = dic[which(dic$PDID==trait),2]
    
    if (length(binwidth_custom)==0){
      binwidth_custom=1
    }
    
    if (Number_of_trials==4){
      trait_name <- df1$DescriptionOfTrait
      
      environ = df1$TRID[1]
      colour <- colour_df[which(colour_df[,1]==environ),2]
      
      p1 <- ggplot(data=df1, aes(x=Score)) + 
        geom_histogram(binwidth = binwidth_custom,
                       fill=colour, 
                       colour= colour, 
                       alpha = 0.8) + 
        labs(title = environ, x=trait_name, y="Count") +
        xlim(c(min-0.5,max+0.5))+
        geom_vline(data=df1, aes(xintercept=median(df1$Score)),
                   linetype="dashed",colour="black",size=0.8) +
        theme_bw()
      
      environ2 = df2$TRID[1]
      colour <- colour_df[which(colour_df[,1]==environ2),2]
      
      p2 <- ggplot(data=df2, aes(x=Score)) + 
        geom_histogram( binwidth = binwidth_custom,
                        fill=colour, 
                        colour= colour, 
                        alpha = 0.8) + 
        labs(title = environ2, x=trait_name, y="Count") +
        xlim(c(min-0.5,max+0.5))+
        geom_vline(data=df1, aes(xintercept=median(df2$Score)),
                   linetype="dashed",colour="black",size=0.8) +
        theme_bw()
      
      environ3 = df3$TRID[1]
      colour <- colour_df[which(colour_df[,1]==environ3),2]
      
      p3 <- ggplot(data=df3, aes(x=Score)) + 
        geom_histogram( binwidth = binwidth_custom,
                        fill=colour, 
                        colour= colour, 
                        alpha = 0.8) + 
        labs(title = environ3, x=trait_name, y="Count") +
        xlim(c(min-0.5,max+0.5))+
        geom_vline(data=df1, aes(xintercept=median(df3$Score)),
                   linetype="dashed",colour="black",size=0.8) +
        theme_bw()
      
      
      environ4 = df4$TRID[1]
      colour <- colour_df[which(colour_df[,1]==environ4),2]
      
      p4 <- ggplot(data=df4, aes(x=Score)) + 
        geom_histogram( binwidth = binwidth_custom,
                        fill=colour, 
                        colour= colour, 
                        alpha = 0.8) + 
        labs(title = environ4, x=trait_name, y="Count") +
        xlim(c(min-0.5,max+0.5))+
        geom_vline(data=df1, aes(xintercept=median(df3$Score)),
                   linetype="dashed",colour="black",size=0.8) +
        theme_bw()
      
      
      
      # combine in one figure and save that figure
      
      PDID = df1$PDID[1]
      filename = paste(path,"PDID",PDID,"_Histograms_",Sys.Date(),".pdf",sep="")
      pdf(filename, height = 7, width =4)
      print(plot_grid(p1, p2, p3,p4, align = "v", nrow = 3))
      
      dev.off()
    }
    


    
    if (Number_of_trials==3){
      trait_name <- df1$DescriptionOfTrait
      
      environ = df1$TRID[1]
      colour <- colour_df[which(colour_df[,1]==environ),2]
      
      p1 <- ggplot(data=df1, aes(x=Score)) + 
        geom_histogram(binwidth = binwidth_custom,
          fill=colour, 
          colour= colour, 
          alpha = 0.8) + 
        labs(title = environ, x=trait_name, y="Count") +
        xlim(c(min-0.5,max+0.5))+
        geom_vline(data=df1, aes(xintercept=median(df1$Score)),
                   linetype="dashed",colour="black",size=0.8) +
        theme_bw()
      
      environ2 = df2$TRID[1]
      colour <- colour_df[which(colour_df[,1]==environ2),2]
      
      p2 <- ggplot(data=df2, aes(x=Score)) + 
        geom_histogram( binwidth = binwidth_custom,
          fill=colour, 
          colour= colour, 
          alpha = 0.8) + 
        labs(title = environ2, x=trait_name, y="Count") +
        xlim(c(min-0.5,max+0.5))+
        geom_vline(data=df1, aes(xintercept=median(df2$Score)),
                   linetype="dashed",colour="black",size=0.8) +
        theme_bw()
      
      environ3 = df3$TRID[1]
      colour <- colour_df[which(colour_df[,1]==environ3),2]
      
      p3 <- ggplot(data=df3, aes(x=Score)) + 
        geom_histogram( binwidth = binwidth_custom,
          fill=colour, 
          colour= colour, 
          alpha = 0.8) + 
        labs(title = environ3, x=trait_name, y="Count") +
        xlim(c(min-0.5,max+0.5))+
        geom_vline(data=df1, aes(xintercept=median(df3$Score)),
                   linetype="dashed",colour="black",size=0.8) +
        theme_bw()
      
      
      # combine in one figure and save that figure
      
      PDID = df1$PDID[1]
      filename = paste(path,"PDID",PDID,"_Histograms_",Sys.Date(),".pdf",sep="")
      pdf(filename, height = 7, width =4)
      print(plot_grid(p1, p2, p3, align = "v", nrow = 3))
      
      dev.off()
     }
    
  
    if (Number_of_trials==2){
      
      trait_name <- df1$DescriptionOfTrait
      
      environ = df1$TRID[1]
      colour <- colour_df[which(colour_df[,1]==environ),2]

      p1 <- ggplot(data=df1, aes(x=Score)) + 
        geom_histogram(binwidth = binwidth_custom,
          fill=colour, 
          colour= colour, 
          alpha = 0.8) + 
        labs(title = environ, x=trait_name, y="Count") +
        xlim(c(min-0.5,max+0.5))+
        geom_vline(data=df1, aes(xintercept=median(df1$Score)),
                   linetype="dashed",colour="black",size=0.8) +
        theme_bw()
      
      environ2 = df2$TRID[1]
      colour <- colour_df[which(colour_df[,1]==environ2),2]
      
      
      p2 <- ggplot(data=df2, aes(x=Score)) + 
        geom_histogram( binwidth = binwidth_custom,
          fill=colour, 
          colour= colour, 
          alpha = 0.8) + 
        labs(title = environ2, x=trait_name, y="Count") +
        xlim(c(min-0.5,max+0.5))+
        geom_vline(data=df1, aes(xintercept=median(df2$Score)),
                   linetype="dashed",colour="black",size=0.8) +
        theme_bw()
      
      # combine in one figure and save that figure
      
      PDID = df1$PDID[1]
      filename = paste(path,"PDID",PDID,"_Histograms_",Sys.Date(),".pdf",sep="")

      pdf(filename, height = 7, width =4)
      print(plot_grid(p1, p2, align = "v", nrow = 3))
      
      dev.off()

    }
    
    if (Number_of_trials==1){
      
      trait_name <- df1$DescriptionOfTrait
      
      environ = df1$TRID[1]
      colour <- colour_df[which(colour_df[,1]==environ),2]
      
      p1 <- ggplot(data=df1, aes(x=Score)) + 
        geom_histogram(binwidth = binwidth_custom,
          fill=colour, 
          colour= colour, 
          alpha = 0.8) + 
        labs(title = environ, x=trait_name, y="Count") +
        xlim(min-0.5,max+0.5)+
        geom_vline(data=df1, aes(xintercept=median(df1$Score)),
                   linetype="dashed",colour="black",size=0.8) +
        theme_bw()
      
      p1
      
      PDID = df1$PDID[1]
      filename = paste(path,"PDID",PDID,"_Histograms_",Sys.Date(),".pdf",sep="")
      
      pdf(filename, height = 7, width =4)
      print(plot_grid(p1, align = "v", nrow = 3))
      
      dev.off()
      
    }
  }
}


# use function ----
histogrammer(df_filt2,"/Users/CathrineKiel/Desktop/test/")





# Outlier removal ----
# When looking at the resulting histograms, outliers are removed manually

df_filt2$Score = as.numeric(as.character(df_filt2$Score))
# plant heigher than 3 m and less than 0.5 m
remove1 <- which(df_filt2$PDID==14 & (df_filt2$Score>300 | df_filt2$Score<50 )) 

# plants that flower for more than 100 days or 0 days 
remove2 <- which((df_filt2$PDID==47 & df_filt2$Score>100) | (df_filt2$PDID==47 & df_filt2$Score==0))

# Seed circularity > 4
remove3 <- which(df_filt2$PDID==23 & df_filt2$Score>4)

# negative date of end flowering
remove4 <- which(df_filt2$PDID==18 & df_filt2$Score<0)

# earliness of flowering >100
remove5 <- which(df_filt2$PDID==17 & df_filt2$Score>100)

# end of flowering that corresponds to duration of 0
ProblemPlots_EndOfFlow <- df_filt2$PLID[which(df_filt2$PDID==47 & df_filt2$Score==0)]
remove5.5 <- which(df_filt2$PDID==18 & df_filt2$PLID %in% ProblemPlots_EndOfFlow)

# maturation date >100
remove6 <- which(df_filt2$PDID==16 & df_filt2$Score>100)

# disease percent >50 as it looks to be a mistake
remove7 <- which(df_filt2$PDID==46 & df_filt2$Score>50)

# disease percent >50 as it looks to be a mistake
remove7.5 <- which(df_filt2$PDID==19 & df_filt2$Score>2000)

# botrytis more than 50 in the late date nordic seed
remove7.6 <- which(df_filt2$PDID==29 & df_filt2$TRID==28.2  & df_filt2$Score>50)





remove_point =c(remove1,remove2,remove4,remove5,remove5.5,remove6,remove7,remove7.5,remove7.6) 

df_filt3 <- df_filt2[-remove_point,] 


# use function again ----
histogrammer(df_filt3,"/Users/CathrineKiel/Desktop/test/")



#Earliness_of_flowering_normalized ----
# now normalize earliness of flowering to dates after sowing 
df_filt11_earliness <- df_filt3[which(df_filt3$PDID==17),]

Scores_ear_date_inJune = df_filt11_earliness$Score-31
Scores_ear_date_inJune_dates <- as.Date(paste("2020","-06-",Scores_ear_date_inJune,sep=""))


for (i in seq(1:length(Scores_ear_date_inJune_dates))){
  
  if (df_filt11_earliness$TRID[i]==29){
    Difference <- Scores_ear_date_inJune_dates[i]-as.Date("2020-04-19") +1
    replacement <-  as.numeric(Difference, units="days")
    df_filt11_earliness$Score[i] = replacement
  }
  
  if (df_filt11_earliness$TRID[i]==28){
    Difference <- Scores_ear_date_inJune_dates[i]-as.Date("2020-04-02") +1
    replacement <-  as.numeric(Difference, units="days")
    df_filt11_earliness$Score[i] = replacement
    }
  if (df_filt11_earliness$TRID[i]==27){
    Difference <- Scores_ear_date_inJune_dates[i]-as.Date("2020-04-09") +1
    replacement <-  as.numeric(Difference, units="days")
    df_filt11_earliness$Score[i] = replacement
  }
  if (df_filt11_earliness$TRID[i]==31){
    Difference <- Scores_ear_date_inJune_dates[i]-as.Date("2020-04-08") +1
    replacement <-  as.numeric(Difference, units="days")
    df_filt11_earliness$Score[i] = replacement
  }
  
}

histogrammer(df_filt11_earliness,"/Users/CathrineKiel/Desktop/test/")








#End_of_flowering_normalized ----
# now normalize end of flowering to dates after sowing 
df_filt11_end <- df_filt3[which(df_filt3$PDID==18),]
df_filt11_end$Date =as.Date(df_filt11_end$Date)
Scores_end_date_inJune = df_filt11_end$Score

for (i in seq(1:length(Scores_end_date_inJune))){
  date = Scores_end_date_inJune[i]
  
  
  if (date <31){
    df_filt11_end$Date[i] =  as.Date(paste("2020","-06-",date,sep=""))
  }
  
  if (date >=31 & date <62){
    date_in_july= date-30
    df_filt11_end$Date[i] =  as.Date(paste("2020","-07-",date_in_july,sep=""))
  }
  
  if (date >=62){
    date_in_aug= date-61
    df_filt11_end$Date[i] =  as.Date(paste("2020","-08-",date_in_aug,sep=""))
    
  }
  
}


for (i in seq(1:nrow(df_filt11_end))){
  
  if (df_filt11_end$TRID[i]==29){
    Difference <- df_filt11_end$Date[i]-as.Date("2020-04-19") +1
    replacement <-  as.numeric(Difference, units="days")
    df_filt11_end$Score[i] = replacement
  }
  
  if (df_filt11_end$TRID[i]==28){
    Difference <- df_filt11_end$Date[i]-as.Date("2020-04-02") +1
    replacement <-  as.numeric(Difference, units="days")
    df_filt11_end$Score[i] = replacement
  }
  if (df_filt11_end$TRID[i]==27){
    Difference <- df_filt11_end$Date[i]-as.Date("2020-04-09") +1
    replacement <-  as.numeric(Difference, units="days")
    df_filt11_end$Score[i] = replacement
  }
  
}

histogrammer(df_filt11_end,"/Users/CathrineKiel/Desktop/test/")







#Maturation_date_normalized ----
# now normalize maturation date to dates after sowing 

df_filt11_mat <- df_filt3[which(df_filt3$PDID==16),]

Scores_end_date_inJuly = df_filt11_mat$Score
df_filt11_mat$Date=as.character(df_filt11_mat$Date)

for (i in seq(1:length(Scores_end_date_inJuly))){
  date = Scores_end_date_inJuly[i]
  
  
  if (date <31){
    df_filt11_mat$Date[i] =  as.Date(paste("2020","-07-",date,sep=""))
  }
  
  if (date >=31 & date <62){
    date_in_july= date-30
    df_filt11_mat$Date[i] =  as.character(paste("2020","-08-",date_in_july,sep=""))
  }
  
  if (date >=62){
    date_in_aug= date-61
    df_filt11_mat$Date[i] =  as.Date(paste("2020","-09-",date_in_aug,sep=""))
    
  }
  
}


for (i in seq(1:nrow(df_filt11_mat))){
  
  if (df_filt11_mat$TRID[i]==29){
    Difference <- as.Date(df_filt11_mat$Date[i])-as.Date("2020-04-19") +1
    replacement <-  as.numeric(Difference, units="days")
    df_filt11_mat$Score[i] = replacement
  }
  
  if (df_filt11_mat$TRID[i]==28){
    Difference <- as.Date(df_filt11_mat$Date[i])-as.Date("2020-04-02") +1
    replacement <-  as.numeric(Difference, units="days")
    df_filt11_mat$Score[i] = replacement
  }
  if (df_filt11_mat$TRID[i]==27){
    Difference <- as.Date(df_filt11_mat$Date[i])-as.Date("2020-04-09") +1
    replacement <-  as.numeric(Difference, units="days")
    df_filt11_mat$Score[i] = replacement
  }
  
}

histogrammer(df_filt11_mat,"/Users/CathrineKiel/Desktop/test/")



# Write resulting file that is quality controlled ----
df_filt11_mat$PDID <- "Maturation_from_sowing"
df_filt11_end$PDID <- "EndOfFlow_from_sowing"
df_filt11_earliness$PDID <- "EarOfFlow_from_sowing"

df_v1 <- rbind(df_filt3,df_filt11_mat,df_filt11_end,df_filt11_earliness)

library(dplyr)
dat <- df_v1 %>%
  mutate(across(everything(), ~ifelse(.=="", NA, as.character(.))))

dat[dat==";"]<-":"


keep = which(colnames(dat) %in% c("PHID","PLID","Rep","Date","Score","Comment_phenotype_Scores","ScoredBy","Comments","PDID","DescriptionOfTrait","Grouping","DescriptionOfMethod","SLID","HarvestDate","HarvestLocation","ParentSLID","GPID","GPSCoordinates","Name","AlternativeName","Donor","GeographicOrigin","Maintaining","TRID","PlotSize","PlotSizeIncludingpaths","SoilType","StartOfTrial","EndOfTrial","Description","Manager","PDID_Name"))
dat_as_desired = dat[,keep]

Minimumformat  = dat[,which(colnames(dat) %in% c("PHID","PLID","Rep","Date","Score","PDID","HarvestDate","HarvestLocation","GPID","Name","AlternativeName","Donor","GeographicOrigin","Maintaining","TRID","StartOfTrial","EndOfTrial","Manager","PDID_Name"))]
write.table(dat_as_desired,paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/","Phenotypes_qualityfiltered",Sys.Date(),".csv",sep="_"),col.names=T, row.names = F, quote=F,sep=",")
write.table(Minimumformat,paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/","Phenotypes_qualityfiltered_minimum",Sys.Date(),".csv",sep="_"),col.names=T, row.names = F, quote=F,sep=",")












# write histogram plots as to include in paper
New_data = Minimumformat
unique(Minimumformat$PDID)

Minimumformat_2 = Minimumformat
listoftraits_renewed = unique(Minimumformat_1$PDID)
listoftraits_renewed_rearranged = listoftraits_renewed[c(9,7,8,10,14,16,1,2,13,11,12,17,19,18,15,3,4,6,5)]

colour_df =as.data.frame(cbind(c(27.0,28.0,28.1,28.2,29.0,31.0),c("#104E8B","#a3172c","#cc9e29","#10108B","#108B8B","pink")))


Minimumformat_2$TRID[which(Minimumformat_2$PDID=="29" & Minimumformat_2$TRID=="28.2")]="28"
Minimumformat_2$TRID[which(Minimumformat_2$PDID=="28" & Minimumformat_2$TRID=="28.1")]="28"


# a function that makes the plots as I desire
for (i in seq(1:length(listoftraits_renewed_rearranged))){
  trait = listoftraits_renewed_rearranged[i]
  
  df2 = Minimumformat_2[which(Minimumformat_2$PDID==trait),]
  
  
  binwidth_custom = dic[which(dic$PDID==trait),2]
  
  if (length(binwidth_custom)==0){
    binwidth_custom=1
  }
  
  # now find plot colours
  trials_present = unique(df2$TRID)
  
  colour <- colour_df[which(colour_df[,1] %in% trials_present),2]
  
  df2$Score = as.numeric(as.character(df2$Score ))
  
  min=NULL
  max=NULL
  
  if (length(min)==0){
    min = min(df2$Score)
  }
  
  if (length(max)==0){
    max = max(df2$Score)
  }
  
  df2$TRID[which(df2$TRID=="27")]="Hor20"
  df2$TRID[which(df2$TRID=="28")]="Dyn20"
  df2$TRID[which(df2$TRID=="28.1")]="Dyn20 (early date)"
  df2$TRID[which(df2$TRID=="28.2")]="Dyn20 (late date)"
  df2$TRID[which(df2$TRID=="29")]="Hor21"
  df2$TRID[which(df2$TRID=="31")]="Dyn21"
  
  df2$TRID = factor(df2$TRID , levels = c("Hor20","Dyn20","Dyn20 (early date)","Dyn20 (late date)","Hor21","Dyn21"))
  
  

  p<- ggplot(data=df2, aes(x=Score,color=TRID,fill=TRID)) + 
    geom_histogram(binwidth = binwidth_custom,alpha=0.3) +
                    #fill=colour, 
                   scale_color_manual(values=colour)  +
                   scale_fill_manual(values=colour)  +
                   #colour= colour, 
                    #alpha = 0.8) + 
    labs(x=trait, y="Count") +
    xlim(c(min-0.5,max+0.5))+
    theme_bw()
  
  p
  
  nam <- paste("p", i, sep = "")
  assign(nam, p)
  
}
Masterplot <- ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19, ncol=4, nrow=5)
Masterplot


Masterplot<-p1 + p2 +p3 +p4 +p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19 + plot_layout(guides = "collect",ncol=3)
save_plot(paste("/Users/CathrineKiel/Desktop/test/",Sys.Date(),".pdf",sep=""), Masterplot,base_width=8,base_height=10,limitsize = FALSE)

