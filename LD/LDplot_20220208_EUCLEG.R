library(ggplot2)
library(tidyverse)
library(devtools)
library(varhandle)

# list of panels
panels=c("EUCLEG")

# list of chr
Chromosomes = c("chr1S","chr1L","chr2","chr3","chr4","chr5","chr6")


# plotting function
LD_plot <- function(panel,chr){
  
    ld_results =NULL
    # assume all
    file.ls2 = intersect(list.files(pattern = ".ld"), list.files(pattern = panel))
    
    for (file in file.ls2){
      
      # if the merged dataset doesn't exist, create it
      if (!exists("ld_results")){
        ld_results <-  read.table(text = gsub(" ", "\t", readLines(file)),head=T,stringsAsFactors = F)
      }
      
      # if the merged dataset does exist, append to it
      if (exists("ld_results")){
        temp_dataset <-read.table(text = gsub(" ", "\t", readLines(file)),head=T,stringsAsFactors = F)
        ld_results<-rbind(ld_results, temp_dataset)
        rm(temp_dataset)
      }
    }
    ld_results$Distance = ld_results$BP_B-ld_results$BP_A
    write.table(ld_results,paste("LDresults",panel,"AllChr",".csv",sep="_"),sep=",",col.names=T,row.names = F,quote=F)
    
    # plot a moving average of LD in increments of 1000,000bp distances
    ld_results_sorted <- ld_results
    ld_results_sorted <- ld_results_sorted[order(ld_results_sorted$Distance),]
    
    bins = seq(1,1000000,by=1000) #5000 bins of 1000 because it is a mapping pop
    
    
    

    # Calculate the average r2 value in each bin.
    my.means=rep(0, length(bins))
    LD.averages=data.frame(bins, my.means,my.means)
    colnames(LD.averages) = c("idx","Avg_distance","Avg_R2")
    
    for (i in 1:length(bins)) {
      data.interval=subset(ld_results, (ld_results$Distance >= bins[i] & ld_results$Distance < (bins[i]+1000000))) 
      
      data.interval = ld_results_sorted[LD.averages$idx[i]:(LD.averages$idx[i]+1000-1),]

      LD.averages$Avg_distance[i]=mean(data.interval$Distance) 
      LD.averages$Avg_R2[i]=mean(data.interval$R2) 
    }
    
    LD.averages=na.omit(LD.averages)
    
    
    # ggplot(LD.averages) +
     #   aes(x = Avg_distance, y = Avg_R2) +
    #   geom_point(shape = "circle", size = 0.5, colour = "#8B8987",alpha=0.01) + 
    #  labs(x = "bp", y = "r2") +
     # theme_classic() +
    #  theme(axis.title.y = element_text(size = 12L, face = "italic"), 
     #       axis.title.x = element_text(size = 12L)) + stat_smooth(aes(outfit=fity <<-..y..,outfit=fitx <<-..x..),method="loess",span=0.1,se=FALSE,color="green",size=0.5) 

    #  maximum <- max(fity) #find maximum LD value on smooth curve
    #  halfmax = maximum/2 # find half of maximum value =LD decay
    # now figure when the smooth curve is equal to halfmax
    # Halfdecay <- round(fitx[min(which(fity<=halfmax))],2) 
    
     Halfdecay<- 681730.4
     halfmax<-0.05512863
    
    filenameplot4=paste("LDplot_bins_smooth",panel,"AllChr",Sys.Date(),".pdf",sep="_")
    p4=ggplot(LD.averages) +
      aes(x = Avg_distance, y = Avg_R2) +
      geom_point(shape = "circle", size = 1.5, colour = "grey50",alpha=0.7) +
      stat_smooth(method = "loess",se=F,colour="green",alpha=1,size=1,span=0.1) +
      #stat_smooth_func(geom="text",method="lm",hjust=0,parse=TRUE) +
      labs(x = "bp (bins of 1000)", y = "r2") +
      geom_hline(yintercept=halfmax, linetype="dashed", color = "black",size=1) +
      xlab(paste("Distance in bp.","LD decay =",Halfdecay,sep=" "))+
      scale_y_continuous(limits = c(0.0,0.25),n.breaks=6) +
      scale_x_continuous(limits = c(0.0,max(LD.averages$Avg_distance)),n.breaks=10) +
      theme_classic()
    
    
    ggsave(filenameplot4, plot = p4, width = 20, height = 14, unit = 'cm')
}
  
# apply function
ld_results=NULL
LD_plot(panel=panels[1],chr=NULL)

