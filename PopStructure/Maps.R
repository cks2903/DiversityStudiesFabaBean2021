library(maps)
#require(devtools)  
#install_github('AndySouth/rworldmap', build_vignettes=TRUE) 
library(rworldmap)
library(RColorBrewer)
library(tidyverse)
library(varhandle)

Mydata = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220316_MapsAndHGW/ForMapping.csv",sep=",",head=T)
head(Mydata)

# subset pop 1
pop1_df = subset(Mydata, Pop == "Pop1")

# count the number in each country
pop1_df$CTR=as.factor(pop1_df$CTR)
Counts_pop1=table(pop1_df$CTR)
Counts_pop1_df=as.data.frame(Counts_pop1)
rm=which(Counts_pop1_df$Var1=="Unknown" | Counts_pop1_df$Var1=="\\N" | Counts_pop1_df$Var1=="0"  | Counts_pop1_df$Freq==0 )
Counts_pop1_df=Counts_pop1_df[-rm,]
colnames(Counts_pop1_df)[1]="country"

visitedMap <- joinCountryData2Map(Counts_pop1_df, 
                                  joinCode = "ISO3",
                                  nameJoinColumn = "country")


pdf("SP1_Continous.pdf", width=20, height=10)
pop1_con <- mapCountryData(visitedMap, 
                       nameColumnToPlot="Freq",
                       catMethod =c(seq(0,25,1)),
                       missingCountryCol ="lightgrey",
                       colourPalette  = colorRampPalette(c("#f0e9f2", "#440154"))(25),
                       addLegend = T,
                       mapTitle = "SP1",
                       border = "black")
dev.off()

# subset pop 2
pop2_df = subset(Mydata, Pop == "Pop2")

# count the number in each country
pop2_df$CTR=as.factor(pop2_df$CTR)
Counts_pop2=table(pop2_df$CTR)
Counts_pop2_df=as.data.frame(Counts_pop2)
rm=which(Counts_pop2_df$Var1=="Unknown" | Counts_pop2_df$Var1=="\\N" | Counts_pop2_df$Var1=="0"  | Counts_pop2_df$Freq==0 )
Counts_pop2_df=Counts_pop2_df[-rm,]
colnames(Counts_pop2_df)[1]="country"

visitedMap <- joinCountryData2Map(Counts_pop2_df, 
                                  joinCode = "ISO3",
                                  nameJoinColumn = "country")



pdf("SP2_Continous.pdf", width=20, height=10)

pop2_con <- mapCountryData(visitedMap, 
                       nameColumnToPlot="Freq",
                       catMethod = c(seq(0,66,1)),
                       missingCountryCol ="lightgrey",
                       colourPalette  = colorRampPalette(c("#fff9c7", "#fde725"))(66),
                       addLegend = T,
                       mapTitle = "SP2",
                       border = "black")
dev.off()

# subset pop 3
pop3_df = subset(Mydata, Pop == "Pop3")

# count the number in each country
pop3_df$CTR=as.factor(pop3_df$CTR)
Counts_pop3=table(pop3_df$CTR)
Counts_pop3_df=as.data.frame(Counts_pop3)
rm=which(Counts_pop3_df$Var1=="Unknown" | Counts_pop3_df$Var1=="\\N" | Counts_pop3_df$Var1=="0"  | Counts_pop3_df$Freq==0 )
Counts_pop3_df=Counts_pop3_df[-rm,]
colnames(Counts_pop3_df)[1]="country"

visitedMap <- joinCountryData2Map(Counts_pop3_df, 
                                  joinCode = "ISO3",
                                  nameJoinColumn = "country")


pdf("SP3_Continous.pdf", width=20, height=10)
pop3_con <- mapCountryData(visitedMap, 
                       nameColumnToPlot="Freq",
                       catMethod = c(seq(0,23,1)),
                       missingCountryCol ="lightgrey",
                       colourPalette  = colorRampPalette(c("#daf2ed", "#21918c"))(23),
                       addLegend = T,
                       mapTitle = "SP3",
                       border = "black")

dev.off()

