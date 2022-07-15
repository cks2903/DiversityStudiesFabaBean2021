# libraries
library(ggplot2)
library(viridis)

# load files
setwd("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220323_FSTplots/")

data1 = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220314_Fst/Full_FST_analyses_weir.fst",head=T,sep="\t")
dim(data1)

# filter

data_filt2=data1



# plots
# note: Pop2 and pop3 appear most differentiated 

for (i in c(4,6,8)){
  current_ana= colnames(data_filt2)[i]
  
  for (chr in c("chr1S","chr1L","chr2","chr3","chr4","chr5","chr6")){
  
    subset=data_filt2[which(data_filt2$CHROM==chr),]
   
    plot=ggplot(subset, aes(x = POS, y = subset[,i])) +
    geom_point(colour="grey",size=3) + ylim(-0.1,1) +
    theme_classic() 
  
    filename=paste("FSTplot",chr,current_ana,Sys.Date(),'.pdf',sep="_")
    ggsave(filename, plot = plot, width = 25, height = 10, unit = 'cm')
  }
}

# try colouring BayeScan markers of interest
SNPs_Sel = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220314_SegregatingMarkersUnderSelection/35_Markers.txt",sep="\t",head=T,stringsAsFactors = F,check.names = F)
highlightlist =colnames(SNPs_Sel)[2:(ncol(SNPs_Sel)-1)]

for (i in c(4,6,8)){
  current_ana= colnames(data_filt2)[i]
  
  subset=data_filt2[which(data_filt2$CHROM=="chr1S"),]
  highlight=subset[which(subset$SNP_name %in% highlightlist),]
  
  if (nrow(highlight)!=0){
      plot=ggplot(subset, aes(x = POS, y = subset[,i])) +
        geom_point(colour="grey",size=3) + ylim(-0.1,1)  +geom_point(data=highlight,
                                                                    aes(x=POS,y=highlight[,i]), colour="green",size=4)+
        theme_classic() 
    
      filename=paste("FSTplot_highlightedBayesScan","chr1S",current_ana,Sys.Date(),'.pdf',sep="_")
      ggsave(filename, plot = plot, width = 25, height = 10, unit = 'cm')
  }
  }

for  (i in c(4,6,8)){
  current_ana= colnames(data_filt2)[i]
  
  subset=data_filt2[which(data_filt2$CHROM=="chr1L"),]
  highlight=subset[which(subset$SNP_name %in% highlightlist),]
  
  if (nrow(highlight)!=0){
    plot=ggplot(subset, aes(x = POS, y = subset[,i])) +
      geom_point(colour="grey",size=3) + ylim(-0.1,1)  +geom_point(data=highlight,
                                                                   aes(x=POS,y=highlight[,i]), colour="green",size=4)+
      theme_classic() 
    
    filename=paste("FSTplot_highlightedBayesScan","chr1L",current_ana,Sys.Date(),'.pdf',sep="_")
    ggsave(filename, plot = plot, width = 25, height = 10, unit = 'cm')
  }
}

for  (i in c(4,6,8)){
  current_ana= colnames(data_filt2)[i]
  
  subset=data_filt2[which(data_filt2$CHROM=="chr2"),]
  highlight=subset[which(subset$SNP_name %in% highlightlist),]
  
  if (nrow(highlight)!=0){
    plot=ggplot(subset, aes(x = POS, y = subset[,i])) +
      geom_point(colour="grey",size=3) + ylim(-0.1,1)  +geom_point(data=highlight,
                                                                   aes(x=POS,y=highlight[,i]), colour="green",size=4)+
      theme_classic() 
    
    filename=paste("FSTplot_highlightedBayesScan","chr2",current_ana,Sys.Date(),'.pdf',sep="_")
    ggsave(filename, plot = plot, width = 25, height = 10, unit = 'cm')
  }
}

for  (i in c(4,6,8)){
  current_ana= colnames(data_filt2)[i]
  
  subset=data_filt2[which(data_filt2$CHROM=="chr3"),]
  highlight=subset[which(subset$SNP_name %in% highlightlist),]
  
  if (nrow(highlight)!=0){
    plot=ggplot(subset, aes(x = POS, y = subset[,i])) +
      geom_point(colour="grey",size=3) + ylim(-0.1,1)  +geom_point(data=highlight,
                                                                   aes(x=POS,y=highlight[,i]), colour="green",size=4)+
      theme_classic() 
    
    filename=paste("FSTplot_highlightedBayesScan","chr3",current_ana,Sys.Date(),'.pdf',sep="_")
    ggsave(filename, plot = plot, width = 25, height = 10, unit = 'cm')
  }
}
for  (i in c(4,6,8)){
  current_ana= colnames(data_filt2)[i]
  
  subset=data_filt2[which(data_filt2$CHROM=="chr4"),]
  highlight=subset[which(subset$SNP_name %in% highlightlist),]
  
  if (nrow(highlight)!=0){
    plot=ggplot(subset, aes(x = POS, y = subset[,i])) +
      geom_point(colour="grey",size=3) + ylim(-0.1,1)  +geom_point(data=highlight,
                                                                   aes(x=POS,y=highlight[,i]), colour="green",size=4)+
      theme_classic() 
    
    filename=paste("FSTplot_highlightedBayesScan","chr4",current_ana,Sys.Date(),'.pdf',sep="_")
    ggsave(filename, plot = plot, width = 25, height = 10, unit = 'cm')
  }
}

for  (i in c(4,6,8)){
  current_ana= colnames(data_filt2)[i]
  
  subset=data_filt2[which(data_filt2$CHROM=="chr5"),]
  highlight=subset[which(subset$SNP_name %in% highlightlist),]
  
  if (nrow(highlight)!=0){
    plot=ggplot(subset, aes(x = POS, y = subset[,i])) +
      geom_point(colour="grey",size=3) + ylim(-0.1,1)  +geom_point(data=highlight,
                                                                   aes(x=POS,y=highlight[,i]), colour="green",size=4)+
      theme_classic() 
    
    filename=paste("FSTplot_highlightedBayesScan","chr5",current_ana,Sys.Date(),'.pdf',sep="_")
    ggsave(filename, plot = plot, width = 25, height = 10, unit = 'cm')
  }
}

for  (i in c(4,6,8)){
  current_ana= colnames(data_filt2)[i]
  
  subset=data_filt2[which(data_filt2$CHROM=="chr6"),]
  highlight=subset[which(subset$SNP_name %in% highlightlist),]
  
  if (nrow(highlight)!=0){
    plot=ggplot(subset, aes(x = POS, y = subset[,i])) +
      geom_point(colour="grey",size=3) + ylim(-0.1,1)  +geom_point(data=highlight,
                                                                   aes(x=POS,y=highlight[,i]), colour="green",size=4)+      theme_classic() 
    
    filename=paste("FSTplot_highlightedBayesScan","chr6",current_ana,Sys.Date(),'.pdf',sep="_")
    ggsave(filename, plot = plot, width = 25, height = 10, unit = 'cm')
  }
}





