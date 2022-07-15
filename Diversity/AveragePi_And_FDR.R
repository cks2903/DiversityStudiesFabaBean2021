library(ggplot2)

file_list <- list.files(pattern="SP1.sites.pi")
Nucleotide_diversity_df <- rep(NA,1000)
Nucleotide_diversity_df=as.data.frame(Nucleotide_diversity_df)
colnames(Nucleotide_diversity_df) = c("Phi_avg")


COUNT=0
for (i in seq(1:length(file_list))){
  df = read.table(file_list[i],head=T)
  Average = mean(df$PI)
  Nucleotide_diversity_df$Phi_avg[i] = Average
  
  if (Average<= 0.26){
    COUNT= COUNT+1
  }

}

FDR= COUNT/1000

filename = paste("Random49_SP1_pi",Sys.Date(),".txt",sep="")
write.table(Nucleotide_diversity_df,filename,quote=F,row.names = F)

Plotname =  paste("Random49_SP1_pi_plot",Sys.Date(),".pdf",sep="")

p = ggplot(data=Nucleotide_diversity_df, aes(x=(as.numeric(Phi_avg)))) +
  geom_histogram(binwidth=0.001,fill = "black",
                 alpha=.8) +
  labs(x="Phi", y=paste("Observations_FDR=",FDR,"_")) +
  geom_vline(xintercept = 0.26,size =0.5, col="green") +
  theme_classic()

ggsave(p,filename=Plotname,height=10,width=15,units="in", dpi=300)




#####################################


file_list <- list.files(pattern="SP2.sites.pi")
Nucleotide_diversity_df <- rep(NA,1000)
Nucleotide_diversity_df=as.data.frame(Nucleotide_diversity_df)
colnames(Nucleotide_diversity_df) = c("Phi_avg")


COUNT=0
for (i in seq(1:length(file_list))){
  df = read.table(file_list[i],head=T)
  Average = mean(df$PI)
  Nucleotide_diversity_df$Phi_avg[i] = Average
  
  if (Average<= 0.26){
    COUNT= COUNT+1
  }
  
}

FDR= COUNT/1000

filename = paste("Random49_SP2_pi",Sys.Date(),".txt",sep="")
write.table(Nucleotide_diversity_df,filename,quote=F,row.names = F)

Plotname =  paste("Random49_SP2_pi_plot",Sys.Date(),".pdf",sep="")

p = ggplot(data=Nucleotide_diversity_df, aes(x=(as.numeric(Phi_avg)))) +
  geom_histogram(binwidth=0.001,fill = "black",
                 alpha=.8) +
  labs(x="Phi", y=paste("Observations_FDR=",FDR,"_")) +
  geom_vline(xintercept = 0.26,size =0.5, col="green") +
  theme_classic()

ggsave(p,filename=Plotname,height=10,width=15,units="in", dpi=300)
