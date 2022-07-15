library(ggplot2)
library(reshape)
library(varhandle)

file = read.table("Genotypes_FarmCPU.txt",sep="\t",head=T,check.names=F)

set.seed(100)
RandomMarkers = sample(colnames(file)[2:ncol(file)],35)
file_filt=file[,c(1,which(colnames(file) %in% RandomMarkers))]

Pop = read.table("pop.txt",sep="\t",head=T,check.names=F)

merged_file = merge(file_filt,Pop,by="Taxa")


meltedfile=melt(merged_file)
meltedfile$variable =unfactor(meltedfile$variable)
meltedfile$Pop =as.factor(meltedfile$Pop)
meltedfile$value =as.factor(meltedfile$value)
meltedfile$variable=factor(meltedfile$variable)


pdf(paste("SegregatingPopsFor35RandomSNPs",Sys.Date(),".pdf",sep=""), 20, 7)
ggplot(meltedfile, aes(Taxa, variable, fill = value, height=0.8)) +
  geom_tile() + scale_fill_manual(values = c("#00C756", "#ff6289", "#213970") )+
  theme_classic()  +facet_grid(~ Pop,scales = "free_x", space = "free")  +  theme(axis.text.x = element_text(angle = 90))
dev.off()
