
library(stringr)
library(ggplot2)
library(RColorBrewer)

# A function that plots multiple GWAS results in one Manhattan plot
# You can specify which files to combine
# and which colours to use for plotting 
# and figure name

MH_multipletrait_plotter <- function(list_of_files, vector_colors, figure_name_with_full_path){
  
  for (i in seq(1:length(list_of_files))){
    
    if (i ==1){
      Results_GWAS = read.table(list_of_files[i],header=T,sep=",")
      trait=str_split(as.character(list_of_files[i]),  "FarmCPU.")[[1]][2]
      trait=str_split(as.character(trait),  ".GWAS")[[1]][1]
      Results_GWAS$trait=trait
    } else{
      
      temp_results =read.table(list_of_files[i],header=T,sep=",")
      trait=str_split(as.character(list_of_files[i]),  "FarmCPU.")[[1]][2]
      trait=str_split(as.character(trait),  ".GWAS")[[1]][1]
      temp_results$trait=trait
      
      Results_GWAS = rbind(Results_GWAS,temp_results,by="SNP")
    }
  }
  
  
  # now plot, 
  Results_GWAS$Chromosome[which(Results_GWAS$Chromosome=="1")]="Chr1S"
  Results_GWAS$Chromosome[which(Results_GWAS$Chromosome=="2")]="Chr1L"
  Results_GWAS$Chromosome[which(Results_GWAS$Chromosome=="3")]="Chr2"
  Results_GWAS$Chromosome[which(Results_GWAS$Chromosome=="4")]="Chr3"
  Results_GWAS$Chromosome[which(Results_GWAS$Chromosome=="5")]="Chr4"
  Results_GWAS$Chromosome[which(Results_GWAS$Chromosome=="6")]="Chr5"
  Results_GWAS$Chromosome[which(Results_GWAS$Chromosome=="7")]="Chr6"
  
  
  Results_GWAS$Position=as.numeric(Results_GWAS$Position)
  Results_GWAS$Chromosome=factor(Results_GWAS$Chromosome, levels= c("Chr1S","Chr1L","Chr2","Chr3","Chr4","Chr5","Chr6"))
  Results_GWAS$P.value=as.numeric(Results_GWAS$P.value)
  Results_GWAS=Results_GWAS[-which(Results_GWAS$trait=="SNP"),]
  Results_GWAS$trait = factor(Results_GWAS$trait,levels=  c(unique(Results_GWAS$trait)))
  
  threshold = -log10(0.05/4790)
  
  

  plot= ggplot()+
    geom_point(data=Results_GWAS, aes(x=Position, y=-log10(P.value), colour = trait),size=4) +
    facet_grid(.~Chromosome,scales="fixed") +
    scale_color_manual(values=vector_colors) +
    geom_hline(yintercept=threshold, linetype="dashed", color = "black") +
    scale_x_continuous(breaks = round(seq(0, 1800000000, by = 150000000),1)) +
    theme(
          axis.text.x = element_text(angle = 90, hjust = 1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")
          )  
  plot
  
#ggsave(figure_name_with_full_path, plot = plot, width = 60, height = 15, unit = 'cm')
ggsave(paste(figure_name_with_full_path,".png",sep=""), plot = plot, width = 60, height = 15, unit = 'cm')

  
}


# Seed trait

list_of_files= c(
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_19.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_27.19.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_29.19.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_31.19.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.TGW_NS_2020_phen_avg.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_20.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_27.20.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_29.20.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_31.20.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.Seed_Area_NS_2020_phen_avg.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_21.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_27.21.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_29.21.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_31.21.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.Seed_Width_NS_2020_phen_avg.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_22.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_27.22.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_29.22.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_31.22.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.Seed_Length_NS_2020_phen_avg.GWAS.Results.csv"
)

vector_colors = c(brewer.pal(n = 6, name = 'Reds')[2:6],brewer.pal(n =6, name = 'Purples')[2:6],brewer.pal(n = 6, name = 'Greens')[2:6],brewer.pal(n = 6, name = 'Blues')[2:6])
figure_name_with_full_path = "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/20220421_OverlayedManhattanplots/Seedtraits"

MH_multipletrait_plotter(list_of_files, vector_colors, figure_name_with_full_path)



# Botrytis
list_of_files= c(
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_29.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_27.29.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_28.2.29.GWAS.Results.csv"
)
vector_colors = c(brewer.pal(n = 5, name = 'Reds')[3],brewer.pal(n =5, name = 'Purples')[3],brewer.pal(n = 5, name = 'Greens')[3])
figure_name_with_full_path = "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/20220421_OverlayedManhattanplots/Botrytis"
MH_multipletrait_plotter(list_of_files, vector_colors, figure_name_with_full_path)


# rust
list_of_files= c(
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_27.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_31.27.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.Rust_.cha.._Se_2021_phen_avg.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_28.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_28.1.28.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.Rust_..._Se_2020_phen_avg.GWAS.Results.csv"
)

vector_colors = c(brewer.pal(n =8, name = 'Reds')[c(3,5,8)],brewer.pal(n =8, name = 'Purples')[c(3,5,8)])
figure_name_with_full_path = "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/20220421_OverlayedManhattanplots/Rust"
MH_multipletrait_plotter(list_of_files, vector_colors, figure_name_with_full_path)


# downy
list_of_files= c(
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_36.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.Downy_Mildew_.cha.._NS_2021_phen_avg.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_29.36.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_46.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.Downy_Mildew_..._NS_2020_early_date_phen_avg.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.Downy_Mildew_..._NS_2020_late_date_phen_avg.GWAS.Results.csv"
)

vector_colors = c(brewer.pal(n =8, name = 'Reds')[c(3,5,8)],brewer.pal(n =8, name = 'Purples')[c(3,5,8)])
figure_name_with_full_path = "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/20220421_OverlayedManhattanplots/Downy"
MH_multipletrait_plotter(list_of_files, vector_colors, figure_name_with_full_path)

# plant height
list_of_files= c(
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_14.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_28.14.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_27.14.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_31.14.GWAS.Results.csv"
)

vector_colors = c(brewer.pal(n = 5, name = 'Reds')[3],brewer.pal(n =5, name = 'Purples')[3],brewer.pal(n = 5, name = 'Greens')[3],brewer.pal(n = 5, name = 'Blues')[3])
figure_name_with_full_path = "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/20220421_OverlayedManhattanplots/Height"
MH_multipletrait_plotter(list_of_files, vector_colors, figure_name_with_full_path)

# Lodging
list_of_files= c(
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_40.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.Lodging_NS_2020_phen_avg.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_31.40.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.Lodging_Se_2021_phen_avg.GWAS.Results.csv"
)

vector_colors = c(brewer.pal(n = 5, name = 'Reds')[3],brewer.pal(n =5, name = 'Purples')[3],brewer.pal(n = 5, name = 'Greens')[3],brewer.pal(n = 5, name = 'Blues')[3])
figure_name_with_full_path = "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/20220421_OverlayedManhattanplots/Lodging"
MH_multipletrait_plotter(list_of_files, vector_colors, figure_name_with_full_path)


# Earliness
list_of_files= c(
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_EarOfFlow_from_sowing.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.EarOfFlow_from_sowing_NS_2020_phen_avg.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.EarOfFlow_from_sowing_Se_2020_phen_avg.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_31.EarOfFlow_from_sowing.GWAS.Results.csv",
  "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/GAPIT.FarmCPU.BLUE_trait_29.EarOfFlow_from_sowing.GWAS.Results.csv"
  
)

vector_colors = c(brewer.pal(n = 5, name = 'Reds')[3],brewer.pal(n =5, name = 'Purples')[3],brewer.pal(n = 5, name = 'Greens')[3],brewer.pal(n = 5, name = 'Blues')[3],brewer.pal(n = 5, name = 'Oranges')[3])
figure_name_with_full_path = "/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/20220421_OverlayedManhattanplots/Earliness"
MH_multipletrait_plotter(list_of_files, vector_colors, figure_name_with_full_path)

