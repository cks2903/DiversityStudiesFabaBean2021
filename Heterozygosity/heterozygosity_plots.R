library(ggplot2)
library(stringr)

Observed_Heterozygosity_means = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/Panel_Characterization/Heterozygosity/ObservedHeterozygosity.txt",sep="\t")
Expected_Heterozygosity = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/Panel_Characterization/Heterozygosity/ExpectedHeterozygosity.txt",sep="\t")

rownames(Observed_Heterozygosity_means)[rownames(Observed_Heterozygosity_means)=="NorFab_MAGIC2"]="MAGIC"
rownames(Observed_Heterozygosity_means)[rownames(Observed_Heterozygosity_means)=="NorFab_Core"]="NorfabCore"
rownames(Expected_Heterozygosity)[rownames(Expected_Heterozygosity)=="NorFab_MAGIC2"]="MAGIC"
rownames(Expected_Heterozygosity)[rownames(Expected_Heterozygosity)=="NorFab_Core"]="NorfabCore"

hetero <-function(panelx){
  
  IndHet <- list.files(path="/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/Panel_Characterization/Heterozygosity/",pattern=paste(panelx,"2022-01-26-1.txt",sep=""))
  hobs_=read.table(IndHet,head=T)
  
  Mean_obs = Observed_Heterozygosity_means$x[which(rownames(Observed_Heterozygosity_means)==panelx)]
  Mean_exp = Expected_Heterozygosity$x[which(rownames(Expected_Heterozygosity)==panelx)]
  colnames(hobs_)=c("Ind","het")
  
  Heterozygosity=ggplot(data=hobs_, aes(x=(as.numeric(as.character(hobs_$het))))) +
    geom_histogram(binwidth=0.01,fill = "black",
                   alpha=.8) +
    labs(x="Multilocus heterozygosity", y="Observations") +
    geom_vline(xintercept = Mean_obs,size =1.5, col="green") +
    geom_vline(xintercept = Mean_exp, size =1.5, col="red")+
    xlim(0,0.4) +theme_classic()
  
  currentDate <- Sys.Date()
  Filename1 <- paste("Heterozygosity",panelx,currentDate,".png",sep="")
  ggsave(Heterozygosity,filename=Filename1,height=10,width=15,units="in", dpi=300)
}

hetero("EUCLEG")
hetero("RSBP")
hetero("KHAZAEI_4PRIL")
hetero("NorfabCore")
hetero("GWB")
hetero("PROFABA")
hetero("MAGIC")
hetero("VICCI")
