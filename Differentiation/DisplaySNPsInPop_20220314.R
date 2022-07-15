library(ggplot2)
library(reshape)
library(varhandle)

file = read.table("5_Markers.txt",sep="\t",head=T)
meltedfile=melt(file)
meltedfile = meltedfile[-1,]

meltedfile$variable =unfactor(meltedfile$variable)
meltedfile$Taxa =unfactor(meltedfile$Taxa)
meltedfile$Pop =as.factor(meltedfile$Pop)
meltedfile$value =as.factor(meltedfile$value)
meltedfile$variable=factor(meltedfile$variable,levels=c("AX.416745027","AX.416737096","AX.416791399","AX.416760427", "AX.416824401"))


pdf(paste("SegregatingPopsFor5CommonSNPs",Sys.Date(),".pdf",sep=""), 20, 7)
ggplot(meltedfile, aes(Taxa, variable, fill = value, height=0.8)) +
  geom_tile() + scale_fill_manual(values = c("#00C756", "#ff6289", "#213970") )+
  theme_classic()  +facet_grid(~ Pop,scales = "free_x", space = "free")  +  theme(axis.text.x = element_text(angle = 90))
dev.off()


file = read.table("35_Markers.txt",sep="\t",head=T,check.names = F)
meltedfile=melt(file)
meltedfile = meltedfile[-1,]

meltedfile$variable =unfactor(meltedfile$variable)
meltedfile$Taxa =unfactor(meltedfile$Taxa)
meltedfile$Pop =as.factor(meltedfile$Pop)
meltedfile$value =as.factor(meltedfile$value)
list = read.delim("order.txt",head=F)
listorder = rev(list[,1])
meltedfile$variable=factor(meltedfile$variable,levels=listorder)


pdf(paste("SegregatingPopsFor35CommonSNPs",Sys.Date(),".pdf",sep=""), 20, 20)
ggplot(meltedfile, aes(Taxa, variable, fill = value, height=0.8)) +
  geom_tile() + scale_fill_manual(values = c("#00C756", "#ff6289", "#213970") )+
  theme_classic()  +facet_grid(~ Pop,scales = "free_x", space = "free")  +  theme(axis.text.x = element_text(angle = 90))
dev.off()

meltedfile <- meltedfile[-which(meltedfile$variable %in% c("AX-416745027","AX-416737096","AX-416791399","AX-416760427", "AX-416824401")),]
pdf(paste("SegregatingPopsFor30CommonSNPsMinus5Best",Sys.Date(),".pdf",sep=""), 20, 20)
ggplot(meltedfile, aes(Taxa, variable, fill = value, height=0.8)) +
  geom_tile() + scale_fill_manual(values = c("#00C756", "#ff6289", "#213970") )+
  theme_classic()  +facet_grid(~ Pop,scales = "free_x", space = "free")  +  theme(axis.text.x = element_text(angle = 90))
dev.off()


                    