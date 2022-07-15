
# Load additional packages not loaded by GAPIT

# GWAS using the FarmCPU implemented in GAPIT

source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")


# Function to do GAPIT GWAS, one at a time
GP_GBLUP<-function(column){
  myGAPIT <- GAPIT(
    Y=myY[,c(1,column)], #fist column is individual ID, the third columns is days to pollination 
    GD=myGD,
    GM=myGM,
    SNP.MAF=0.05,
    PCA.total=3,
    model=c("FarmCPU")
    )
}
print("GAPIT function loaded")


# Import data
myY <- read.table("../Pheno.csv", head = TRUE, sep=",") # Load phenotypes
print("Phenotypes loaded")
myGD <- read.table("../Genotypes_FarmCPU.txt", head = TRUE) # Load genotypes, missingness is not allowed
print("Genotypes loaded")
myGM <- read.table("../GenotypeInfo_FarmCPU.txt", head = TRUE)
print("Genotype information loaded")



for (i in (2:ncol(myY))){
  GP_GBLUP(i)
}
