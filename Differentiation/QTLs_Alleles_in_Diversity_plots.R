
library(ggplot2)
library(varhandle)
library(cowplot)
library(data.table)
library(ggpubr)
library(rstatix)

setwd("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20221104_FavorableAlles")

Info = read.table("LowPhenotypeAlleles.csv",sep=",",head=T)
Info = Info[1:65,]
colnames(Info)[2]="ID"

# load original genotypes
og_genotypes <- fread("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/Raw_files/MakingFinal_vcfFile_20220106/Mapped_HighQualSNPs2022-01-06.vcf")
og_genotypes_df=as.data.frame(og_genotypes) 

# diversity panel with their pops
popinfo = read.table("Popinfo.csv",sep=",",head=T)


# now we calculate the frequency of the favorable allele in each pop
og_genotypes_df_filt=og_genotypes_df[c(1:9,which(colnames(og_genotypes_df) %in% popinfo[,3]))]

merged = merge(og_genotypes_df_filt,Info, all.y=T,by="ID")

merged_t = t(merged)

Dataframe = merged_t
colnames(Dataframe) = Dataframe[1,]
Dataframe=Dataframe[10:nrow(Dataframe),]
Dataframe=as.data.frame(Dataframe)
Dataframe$CELcode =rownames(Dataframe)

merged2 = merge(Dataframe,popinfo,by="CELcode")

for (i in seq(2,68)){
  
  #Mean=table(merged2[,i])
  #data_for_SNP = cbind(table(merged2[,i],merged2$ADMIXTURE.Subpopulation),Mean)
  data_for_SNP = table(merged2[,i],merged2$ADMIXTURE.Subpopulation)
  
  if (length(which(rownames(data_for_SNP)=="./."))!=0){
    data_for_SNP=data_for_SNP[-which(rownames(data_for_SNP)=="./."),]
  }
  SNPname = colnames(merged2)[i]
  AlleleWithNefgativeEffect = Info$Allele.with.negative.effect.on.genotype[which(Info$ID==SNPname)]
  Trait = Info$Trait[which(Info$ID==SNPname)]
  
  
  data_for_SNP_=melt(data_for_SNP)
  #print(paste(i,data_for_SNP_$Var1[2]))
  data_for_SNP_$Var1=unfactor(data_for_SNP_$Var1)
  
  if (nrow(data_for_SNP_)==9){
    
    data_for_SNP_$value[1]=data_for_SNP_$value[1]+0.5*data_for_SNP_$value[2]
    data_for_SNP_$value[3]=data_for_SNP_$value[3]+0.5*data_for_SNP_$value[2]
    
    data_for_SNP_$value[4]=data_for_SNP_$value[4]+0.5*data_for_SNP_$value[5]
    data_for_SNP_$value[6]=data_for_SNP_$value[6]+0.5*data_for_SNP_$value[5]
    
    data_for_SNP_$value[7]=data_for_SNP_$value[7]+0.5*data_for_SNP_$value[8]
    data_for_SNP_$value[9]=data_for_SNP_$value[9]+0.5*data_for_SNP_$value[8]
    
    data_for_SNP_=data_for_SNP_[-c(2,5,8),]
  }
  
  data_for_SNP_$Var1[which(data_for_SNP_$Var1==paste(AlleleWithNefgativeEffect,"/",AlleleWithNefgativeEffect,sep=""))]="Negative Allele"

  data_for_SNP_$Var1[which(data_for_SNP_$Var1!="Negative Allele")]="Positive Allele"
  
  
  # statistics
  
  #contingency table
  data_for_SNP
  data_for_SNP_
  Reshaped = reshape(data_for_SNP_, idvar = "Var2", timevar = "Var1", direction = "wide")
  
  SP1vsSP2_Contingency = Reshaped[c(which(Reshaped$Var2=="SP1"),which(Reshaped$Var2=="SP2")),]
  SP1vsSP3_Contingency = Reshaped[c(which(Reshaped$Var2=="SP1"),which(Reshaped$Var2=="SP3")),]
  SP2vsSP3_Contingency = Reshaped[c(which(Reshaped$Var2=="SP2"),which(Reshaped$Var2=="SP3")),]
  
  Threeasterisks = 0.005/(65*3)
  Twoasterisks = 0.01/(65*3)
  Oneasterisks = 0.05/(65*3)
  
  testSP1vsSP2 = fisher.test(data.matrix(SP1vsSP2_Contingency))
  pvalue_testSP1vsSP2 = testSP1vsSP2$p.value
  
  testSP1vsSP3 = fisher.test(data.matrix(SP1vsSP3_Contingency))
  pvalue_testSP1vsSP3 = testSP1vsSP3$p.value
  
  testSP2vsSP3 = fisher.test(data.matrix(SP2vsSP3_Contingency))
  pvalue_testSP2vsSP3 = testSP2vsSP3$p.value
  
  if (pvalue_testSP1vsSP2<=Threeasterisks){
    pvalue_testSP1vsSP2 ="***"
  }
  
  if (pvalue_testSP1vsSP2<=Twoasterisks & pvalue_testSP1vsSP2>=Threeasterisks){
    pvalue_testSP1vsSP2 ="**"
  }
  
  if (pvalue_testSP1vsSP2<=Oneasterisks & pvalue_testSP1vsSP2>=Twoasterisks){
    pvalue_testSP1vsSP2 ="*"
  }
  
  if (pvalue_testSP1vsSP2>=Oneasterisks){
    pvalue_testSP1vsSP2 ="NS"
  }
  
  if (pvalue_testSP1vsSP3<=Threeasterisks){
    pvalue_testSP1vsSP3 ="***"
  }
  
  if (pvalue_testSP1vsSP3<=Twoasterisks & pvalue_testSP1vsSP3>=Threeasterisks){
    pvalue_testSP1vsSP3 ="**"
  }
  
  if (pvalue_testSP1vsSP3<=Oneasterisks & pvalue_testSP1vsSP3>=Twoasterisks){
    pvalue_testSP1vsSP3 ="*"
  }
  
  if (pvalue_testSP1vsSP3>=Oneasterisks){
    pvalue_testSP1vsSP3 ="NS"
  }
  
  if (pvalue_testSP2vsSP3<=Threeasterisks){
    pvalue_testSP2vsSP3 ="***"
  }
  
  if (pvalue_testSP2vsSP3<=Twoasterisks & pvalue_testSP2vsSP3>=Threeasterisks){
    pvalue_testSP2vsSP3 ="**"
  }
  
  if (pvalue_testSP2vsSP3<=Oneasterisks & pvalue_testSP2vsSP3>=Twoasterisks){
    pvalue_testSP2vsSP3 ="*"
  }
  
  if (pvalue_testSP2vsSP3>=Oneasterisks){
    pvalue_testSP2vsSP3 ="NS"
  }
  
  
  
  #plotting
  
 plot=ggplot(data_for_SNP_, aes(x =Var2, y=value, fill =Var1)) + ylab("Allele frequency")+
  geom_bar(position="fill", stat="identity") +
  
   annotate("text", x = 1.5, y = 1.2, 
            label = pvalue_testSP1vsSP2,
            size = 12, color = "black") +
   
   annotate("text", x = 2, y = 1.4, 
            label = pvalue_testSP1vsSP3,
            size = 12, color = "black") +
   
   annotate("text", x = 2.5, y = 1.2, 
            label = pvalue_testSP2vsSP3,
            size = 12, color = "black") +
   theme_classic() +xlab(colnames(merged2)[i])+ ggtitle(Trait) +scale_fill_manual(values=c("pink","navy")) + theme(legend.position="none",text = element_text(size = 40))
 
 nm <- paste("p",i,sep="") # name of your variable
 assign(nm,plot)
 
ggsave(paste(nm,".png",sep=""),plot,width=18,height=6)
        
}



#megaplot=plot_grid(p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,
     #     p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,
      #    p31,p32,p33,p34,p35,p36,p37,p38,p39,
       #   p40,p41,p42,p43,p44,p45,p46,p47,
        #  p48,p49,p50,p51,p52,p53,p54,p55,p56,p57,p58,p59,p60,
         # p61,p62,p63,p64,p65,nrow=13,align="v")




