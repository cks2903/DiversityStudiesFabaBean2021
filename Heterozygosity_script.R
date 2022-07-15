##################################################################################
############ Calculation of the inbreeding coefficients of core lines ############ 
##################################################################################

# load libraries needed
#install.packages("inbreedR")
library("inbreedR")
library("ggplot2")
library("agricolae")
library("adegenet")
library("hierfstat")
library("pegas")
library("varhandle")
args=commandArgs(trailingOnly = TRUE)


# load functions
{
  is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
}

# load data
{
genotypes=read.table("/home/cks/norfab/faststorage/cks/20220120_DiversityPanel/20220120_Panel_Specific_Analyses/20220126_RedundancyWithinPanels/Mapped_HighQualSNPs2022_RedFiltered20220126.vcf", header = T, sep = "\t",stringsAsFactors = F) #vcf file, not MAF filtered
dim(genotypes)
genotypes_shorter=genotypes[,10:ncol(genotypes)]
genotypes_t=t(genotypes_shorter)
dim(genotypes_t)
rownames(genotypes_t)=colnames(genotypes_shorter)
colnames(genotypes_t)=genotypes[,3]
genotypes=genotypes_t
names=rownames(genotypes)
}

# format data 
{
genotypes_conv=genotypes
genotypes_conv[genotypes_conv=="1|0"]=1
genotypes_conv[genotypes_conv=="0|1"]=1
genotypes_conv[genotypes_conv=="0|0"]=0
genotypes_conv[genotypes_conv=="1|1"]=0 #Changing genotype format, so all heterozygotes are 1, and homozygotes are 0

# check if data has the right format
num_ind=as.numeric(as.character(args[2]))
num_loci=as.numeric(as.character(args[3]))

dim(genotypes_conv)
genotypes_conv[1:10,1:10]
genotypes_conv=as.data.frame(genotypes_conv) 


for (i in seq(1:ncol(genotypes_conv))){
  genotypes_conv[,i] = unfactor(genotypes_conv[,i])
  genotypes_conv[,i] = as.numeric(genotypes_conv[,i])
}
dim(genotypes_conv)

str(genotypes_conv)
#check_data(genotypes_conv, num_ind = num_ind, num_loci = num_loci)
genotypes_conv=as.data.frame(genotypes_conv)
rownames(genotypes_conv)=names
}


# divide df for each panel
Dic <- read.table("/home/cks/norfab/faststorage/cks/20220120_DiversityPanel/20220120_Panel_Specific_Analyses/20220126_RedundancyWithinPanels/Lines_Panel_DicNew.txt",head=T,sep="\t")


Panels = unique(Dic[,2])

Accessions_EUCLEG = Dic[which(Dic[,2]=="EUCLEG"),1]
Accessions_RSBP = Dic[which(Dic[,2]=="RSBP"),1]
Accessions_Fourwaycross = Dic[which(Dic[,2]=="Four-way-cross"),1]
Accessions_NORFAB = Dic[which(Dic[,2]=="NORFAB"),1]
Accessions_GWB = Dic[which(Dic[,2]=="GWB"),1]
Accessions_ProFaba = Dic[which(Dic[,2]=="ProFaba"),1]
Accessions_SevenparentMAGIC = Dic[which(Dic[,2]=="Seven-parent-Seven-parent-MAGIC"),1]
Accessions_VICCI = Dic[which(Dic[,2]=="VICCI"),1]


df_EUCLEG <- genotypes_conv[which(rownames(genotypes_conv) %in% Accessions_EUCLEG),]
df_RSBP <- genotypes_conv[which(rownames(genotypes_conv) %in% Accessions_RSBP),]
df_Fourwaycross <- genotypes_conv[which(rownames(genotypes_conv) %in% Accessions_Fourwaycross),]
df_NORFAB <- genotypes_conv[which(rownames(genotypes_conv) %in% Accessions_NORFAB),]
df_GWB <- genotypes_conv[which(rownames(genotypes_conv) %in% Accessions_GWB),]
df_ProFaba<- genotypes_conv[which(rownames(genotypes_conv) %in% Accessions_ProFaba),]
df_SevenparentMAGIC<- genotypes_conv[which(rownames(genotypes_conv) %in% Accessions_SevenparentMAGIC),]
df_VICCI<- genotypes_conv[which(rownames(genotypes_conv) %in% Accessions_VICCI),]



# Plot one panel at a time

hetero <-function(data,panelx){

	het <- MLH(data)
	Heterozygositylevels=cbind(rownames(data),het)
	Heterozygositylevels=as.data.frame(Heterozygositylevels)

	Heterozygosity=ggplot(data=Heterozygositylevels, aes(x=(as.numeric(as.character(Heterozygositylevels$het))))) +
  geom_histogram(binwidth=0.01,fill = "black",
                 alpha=.8) +
  labs(x="Multilocus heterozygosity", y="Observations") +
  geom_vline(xintercept = mean(as.numeric(as.character(Heterozygositylevels$het,na.rm=T))),size =1.5, col="green") +
  xlim(0,0.4) +theme_classic()

  currentDate <- Sys.Date()
	Filename1 <- paste("Heterozygosity",panelx,currentDate,".png",sep="")
	ggsave(Heterozygosity,filename=Filename1,height=10,width=15,units="in", dpi=300)


Filename2 <- paste("Heterozygosity",panelx,currentDate,".txt",sep="")
write.table(Heterozygositylevels,Filename2,row.names=F, col.names=T, sep="\t",quote=F) 

}

hetero(df_EUCLEG,"EUCLEG")
hetero(df_RSBP,"RSBP")
hetero(df_Fourwaycross,"Four-way-cross")
hetero(df_NORFAB,"NORFAB")
hetero(df_GWB,"GWB")
hetero(df_ProFaba,"ProFaba")
hetero(df_SevenparentMAGIC,"Seven-parent-MAGIC")
hetero(df_VICCI,"df_VICCI")

