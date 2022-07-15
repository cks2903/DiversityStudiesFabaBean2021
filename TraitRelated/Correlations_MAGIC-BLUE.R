library(corrgram)
library(corrplot)
library(agricolae)

df = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220527_MAGIC_BLUP_Cors/20220530_BLUEs_for_cor.csv",sep=",",head=T,colClasses = c("character",rep("numeric",55))) #File with phenotypes used for GWAS
dim(df)

df_filtl =df[,2:ncol(df)]
correlations <- cor(df_filtl,use = "na.or.complete")

p1= corrgram(df_filtl,
         main="Trait correlations",
         lower.panel=panel.pts, upper.panel=panel.conf,
         diag.panel=panel.density)


p2= corrgram(df_filtl,
             main="Trait correlations",
             upper.panel=panel.cor,
             diag.panel=panel.density)
p2

filename1= paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220527_MAGIC_BLUP_Cors/","correlations_all_traits_reordered",Sys.Date(),".pdf",sep="")
pdf(file = filename1)
p3 = corrplot(cor(df_filtl,use = "na.or.complete"), addCoef.col = 1,number.cex = 0.2, tl.cex = 0.5)    # Draw corrplot with default font size
p3
dev.off()

pval_threshold = 0.05/((55*55-55)/2)
res1 <- cor.mtest(df_filtl, conf.level = .95)
filename15= paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220527_MAGIC_BLUP_Cors/","correlations_all_traits_reordered_sig",Sys.Date(),".pdf",sep="")
pdf(file = filename15)
p15 = corrplot(cor(df_filtl,use = "na.or.complete"),sig.level = pval_threshold,pch.cex = 0.5,pch.col = "green",insig = "label_sig", addCoef.col = 1,number.cex = 0.2, tl.cex = 0.5,p.mat=res1$p)    # Draw corrplot with default font size
p15
dev.off()


# vælg ud hvor LD group1 er significant
keep = which(colnames(df_filtl) %in% c("Lodging..Hor21.","Susc..to.dm.perc...Dyn20..late.date.","Susc..to.dm.perc...All.trials.","Seed.length..Hor20.","Seed.length..All.trials.","Seed.area..Dyn21.","Seed.width..Hor21." ,"Seed.width..All.trials.","Seed.width..Hor20.","Seed.width..Dyn21.","TGW..Hor20.","Seed.area..All.trials."  ,"Seed.area..Dyn20.","Seed.width..Dyn20.","Seed.area..Hor20.","Susc..to.rust.perc...All.trials.","Plant.height..Dyn21.","Susc..to.rust.perc...Dyn20..early.date.","End.of.flowering..Dyn20.","Plant.height..All.trials.","Plant.height..Hor20.","Susc..to.cs..All.trials.","Sterile.tillers..Dyn20.","Susc..to.cs..Hor20." ))

df_filt2 = df_filtl[,keep]
dim(df_filt2)


pval_threshold = 0.05/((24*24-24)/2)
res1 <- cor.mtest(df_filt2, conf.level = .95)
filename25= paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220527_MAGIC_BLUE_Cors/","correlations_LDGroup1_reordered_significant",Sys.Date(),".pdf",sep="")
pdf(file = filename25)
p25 = corrplot(cor(df_filt2,use = "na.or.complete"),sig.level = pval_threshold,pch.cex = 1,pch.col = "green",insig = "label_sig", addCoef.col = 1,number.cex = 0.35, tl.cex = 1,p.mat=res1$p)    # Draw corrplot with default font size
p25
dev.off()




filename2= paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220527_MAGIC_BLUE_Cors/","correlations_LDGroup1_reordered",Sys.Date(),".pdf",sep="")
pdf(file = filename2)
p4 = corrplot(cor(df_filt2,use = "na.or.complete"), addCoef.col = 1,number.cex =0.35, tl.cex = 0.9)    # Draw corrplot with default font size
p4
dev.off()


cors = correlation(df_filtl)
cors1 = cbind(cors$correlation,cors$pvalue)




