library(tidyverse)
pca <- read_table2("pca.eigenvec", col_names = FALSE)
eigenval <- scan("pca.eigenval")

# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))


# get panel information
Dic <- read.table("/home/cks/norfab/faststorage/cks/20220120_DiversityPanel/20220120_Panel_Specific_Analyses/20220120_LD/Lines_Panel_DicNew.txt",head=F,sep="\t")
colnames(Dic)[1]="ind"

pca_copy = as.data.frame(pca)
merged <- merge(pca_copy,Dic,by="ind",sort = FALSE)
all(merged[,1]==pca[,1])

#remake
pca <- as.tibble(data.frame(pca, merged$V2))

#plotting

# first convert to percentage variance explained
pve <- data.frame(PC = 1:nrow(pca), pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
b <- a + ylab("Percentage variance explained") + theme_light()

ggsave("AllPanels_PercVarExplained.png", plot = b, width = 20, height = 14, unit = 'cm')



# plot pca
min= min(c(min(pca$PC1),min(pca$PC2)))
max =max(c(max(pca$PC1),max(pca$PC2)))


c <- ggplot(pca, aes(PC1, PC2, col =merged.V2 )) + geom_point(size = 1)
c <- c +   scale_colour_manual(values = c("#0D2C54","#00AF54","#C6878F")) 
c <- c + coord_equal() + theme_classic() + xlim(min,max) +ylim(min,max) 
d<-c + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))
ggsave("PC1PC2All.png", plot = d, width = 20, height = 15, unit = 'cm')


min= min(c(min(pca$PC2),min(pca$PC3)))
max =max(c(max(pca$PC2),max(pca$PC3)))

e <- ggplot(pca, aes(PC2, PC3, col =merged.V2 )) + geom_point(size = 1)
e <- e +   scale_colour_manual(values = c("#0D2C54","#00AF54","#C6878F")) 
e <- e + coord_equal() + theme_classic()  + xlim(min,max) +ylim(min,max) 
f<-e + xlab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) + ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))
ggsave("PC2PC3All.png", plot = f, width = 20, height = 15, unit = 'cm')

