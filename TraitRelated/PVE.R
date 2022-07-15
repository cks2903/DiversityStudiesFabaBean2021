# libraries
library("pheatmap")
library(RColorBrewer)

# Load observations, all traits
Phenotypes = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/_Phenotypes_qualityfiltered_minimum_2022-03-25_.csv",sep=",",head=T) #load phenotype data
Phenotypes$PDIDTRID = paste(Phenotypes$PDID,Phenotypes$TRID,sep=":")

# Load genotypes
Genotypes1 = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/GWAS/Genotypes_FarmCPU.txt",sep="\t",head=T, check.names = F) #load genotypes, filtered at MAF 5%
colnames(Genotypes1)

# Load list of the 35 selection markers
SelMarkers = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220419_35RandomMarkersStudy/35MarkersUndersel.txt",sep="\t",head=F, check.names=F) # dataframe with 35 selection markers, col1 includes SNP name, col2 chromosome and col3 physical position
head(SelMarkers)
SNPs = SelMarkers$V1

LDBlock1 = c("AX-416824401","AX-416760427","AX-416791399","AX-416783057","AX-181165197","AX-416724016","AX-181497981")
LDBlock2 = c("AX-181482613","AX-416771656","AX-181188041","AX-181487950","AX-181486832")
LDBlock3 = c("AX-181487950","AX-416767699")

# GWAS dictionary
GWASdic= read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/Trait_GWAS_fileNames.txt",sep="\t",head=T) # a directory that associate each trait with the name of the FarmCPU produced GWAS results file

# For each of the 35 markers, the variance explained of all traits should be tested. 


All_traits = c(unique(Phenotypes$PDID),unique(Phenotypes$PDIDTRID))


PVE_function <- function(trait,GWAS=0){
  
  PercentageVarianceExplained_20RandomMarkersSUM = 0
  FDR_COUNT_20SelectionMarkers = 0
  FDR_COUNT_20GWASMarkers = 0
  PercentageVarianceExplained_1RandomMarkersSUM = 0
  FDR_COUNT_5SelectionMarkers = 0
  FDR_COUNT_3SelectionMarkers = 0
  FDR_COUNT_2SelectionMarkers = 0
  PercentageVarianceExplained_5RandomMarkersSUM = 0
  PercentageVarianceExplained_2RandomMarkersSUM = 0
  PercentageVarianceExplained_3RandomMarkersSUM = 0
  
  
  

  Results_df<- data.frame(matrix(ncol = 51, nrow = 1))
  colnames(Results_df) = c("trait",SNPs)
  Results_df[1]<- trait
  
  Random_marker_df<- data.frame(matrix(ncol = 1001, nrow = 1))
  colnames(Random_marker_df) = c("trait",seq(1:1000))
  Random_marker_df[1]<- trait
  
  Phenotype_relevant <- Phenotypes[which(Phenotypes$PDID==trait | Phenotypes$PDIDTRID==trait),]
  Phenotype_relevant =Phenotype_relevant[,c(8,3)] # lines and score
  colnames(Phenotype_relevant)=c("Name","Score")
  
  
  #How much can one random marker explain?
  
  for (p in seq(1:1000)){
    Random_Marker  = sample(colnames(Genotypes1)[2:ncol(Genotypes1)],1)
    Genotype_relevant <- Genotypes1[,c(1,which(colnames(Genotypes1)==Random_Marker))]
    colnames(Genotype_relevant)[1]="Name"
    merged = merge(Phenotype_relevant,Genotype_relevant,by="Name",all = T)
    merged_filt = na.omit(merged)
    colnames(merged_filt)[1]=c("Line")
    merged_filt$Line=as.factor(merged_filt$Line)
    merged_filt$Score=as.numeric(as.character((merged_filt$Score)))
  
    colnames(merged_filt)[3:ncol(merged_filt)]=paste("SNP",seq(1:20),sep="")
    # Check variance explained by SNP
  
    model = lm(Score ~SNP1 ,data=merged_filt)
    PredictedValues = predict(model)
  
    Actualmean = mean(na.omit(merged_filt$Score))
    Predictedmean= mean(PredictedValues)
  
    top = 0
    bottom = 0
  
    for (i in 1:nrow(merged_filt)){
    
      add_top = (PredictedValues[i]-Predictedmean)^2
      add_bottom = (merged_filt$Score[i]-Actualmean)^2
    
      top = top + add_top
      bottom = bottom + add_bottom
    }
  
    r2 = top /bottom
    PercentageVarianceExplained_OneRandomMarker= round(r2*100,1)
  
    Random_marker_df[p+1]=PercentageVarianceExplained_OneRandomMarker
    Results_df[41] = PercentageVarianceExplained_OneRandomMarker
  }
  
  
  
  
  
  # PVE explained for each marker under selection
  for (k in seq(1:35)){
    SNP=SNPs[k]
    Genotype_relevant <- Genotypes1[,c(1,which(colnames(Genotypes1)==SNP))]
    
    if (length(which(colnames(Genotypes1)==SNP))==0){
      next
    }
  
    colnames(Genotype_relevant)[1]="Name"
    merged = merge(Phenotype_relevant,Genotype_relevant,by="Name",all = T)
    merged_filt = na.omit(merged)
    colnames(merged_filt)=c("Line","Score","Marker")
    #merged_filt$Marker=as.factor(merged_filt$Marker)
    merged_filt$Line=as.factor(merged_filt$Line)
    merged_filt$Score=as.numeric(as.character((merged_filt$Score)))
    
    # Check variance explained by SNP
    model = lm(Score ~Marker ,data=merged_filt)
    PredictedValues = predict(model)
    
    Actualmean = mean(na.omit(merged_filt$Score))
    Predictedmean= mean(PredictedValues)
    
    top = 0
    bottom = 0
    
    for (i in 1:nrow(merged_filt)){
      
      add_top = (PredictedValues[i]-Predictedmean)^2
      add_bottom = (merged_filt$Score[i]-Actualmean)^2
      
      top = top + add_top
      bottom = bottom + add_bottom
    }
    
    r2 = top /bottom
    PercentageVarianceExplained = round(r2*100,1)
    
    Results_df[(k+1)] <- PercentageVarianceExplained
    
    
    
  }
  
  #Then see how much all selection markers can explain
  Genotype_relevant <- Genotypes1[,c(1,which(colnames(Genotypes1)%in% SNPs))]
  colnames(Genotype_relevant)[1]="Name"
  merged = merge(Phenotype_relevant,Genotype_relevant,by="Name",all = T)
  merged_filt = na.omit(merged)
  colnames(merged_filt)[1]=c("Line")
  merged_filt$Line=as.factor(merged_filt$Line)
  merged_filt$Score=as.numeric(as.character((merged_filt$Score)))
  
  colnames(merged_filt)[3:ncol(merged_filt)]=paste("SNP",seq(1:20),sep="")
  # Check variance explained by SNP
  
  model = lm(Score ~SNP1+SNP2+SNP3+SNP4+SNP5+SNP6+SNP7+SNP8+SNP9+SNP10+SNP11+SNP12+SNP13+SNP14+SNP15+SNP16+SNP17+SNP18+SNP19+SNP20 ,data=merged_filt)

  PredictedValues = predict(model)
  
  Actualmean = mean(na.omit(merged_filt$Score))
  Predictedmean= mean(PredictedValues)
  
  top = 0
  bottom = 0
  
  for (i in 1:nrow(merged_filt)){
    
    add_top = (PredictedValues[i]-Predictedmean)^2
    add_bottom = (merged_filt$Score[i]-Actualmean)^2
    
    top = top + add_top
    bottom = bottom + add_bottom
  }
  
  r2 = top /bottom
  PercentageVarianceExplained_AllMarkers = round(r2*100,1)
  
  Results_df[37] =PercentageVarianceExplained_AllMarkers
  
  
  # check PVE explained by top 20 GWAS markers, only run if argument GWAS is set to 1, by default 0
  if (GWAS==1){
    GWASfilename = GWASdic[which(GWASdic$trait==trait),2]
    GWASfilenamefullpath = paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/",GWASfilename,sep="")
    files_to_look_in = list.files("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/Results/GWAS_20220328/PC3/")
    if(GWASfilename %in% files_to_look_in){
      GWAS_results = read.table(GWASfilenamefullpath,header=T,sep=",")
      
      
      TopMarkers = GWAS_results$SNP[1:20]
      
      Genotype_relevant <- Genotypes1[,c(1,which(colnames(Genotypes1)%in% TopMarkers))]
      colnames(Genotype_relevant)[1]="Name"
      merged = merge(Phenotype_relevant,Genotype_relevant,by="Name",all = T)
      merged_filt = na.omit(merged)
      colnames(merged_filt)[1]=c("Line")
      merged_filt$Line=as.factor(merged_filt$Line)
      merged_filt$Score=as.numeric(as.character((merged_filt$Score)))
      
      colnames(merged_filt)[3:ncol(merged_filt)]=paste("SNP",seq(1:20),sep="")
      # Check variance explained by all 20 random SNPs
      
      model1 = lm(Score ~SNP1+SNP2+SNP3+SNP4+SNP5+SNP6+SNP7+SNP8+SNP9+SNP10+SNP11+SNP12+SNP13+SNP14+SNP15+SNP16+SNP17+SNP18+SNP19+SNP20,data=merged_filt)
      PredictedValues = predict(model1)
      
      Actualmean = mean(na.omit(merged_filt$Score))
      Predictedmean= mean(PredictedValues)
      
      top = 0
      bottom = 0
      
      for (i in 1:nrow(merged_filt)){
        
        add_top = (PredictedValues[i]-Predictedmean)^2
        add_bottom = (merged_filt$Score[i]-Actualmean)^2
        
        top = top + add_top
        bottom = bottom + add_bottom
      }
      
      
      r2 = top /bottom
      PercentageVarianceExplained_20TopGWASMarkers = round(r2*100,1)
      Results_df[40] =PercentageVarianceExplained_20TopGWASMarkers
      
    }
  }
  
  #Then see how 20 random non-selection markers can explain, repeat 1000 times
  for (h in seq(1:1000)){
    
    SNPs_not_under_Selection <- colnames(Genotypes1)[-which(colnames(Genotypes1)%in% SNPs)]
    SNPs_not_under_Selection = SNPs_not_under_Selection[2:length(SNPs_not_under_Selection)]
    random_marker =sample(SNPs_not_under_Selection,20)
    Genotype_relevant <- Genotypes1[,c(1,which(colnames(Genotypes1)%in% random_marker))]
    colnames(Genotype_relevant)[1]="Name"
    merged = merge(Phenotype_relevant,Genotype_relevant,by="Name",all = T)
    merged_filt = na.omit(merged)
    colnames(merged_filt)[1]=c("Line")
    merged_filt$Line=as.factor(merged_filt$Line)
    merged_filt$Score=as.numeric(as.character((merged_filt$Score)))
  
    colnames(merged_filt)[3:ncol(merged_filt)]=paste("SNP",seq(1:20),sep="")
    # Check variance explained by all 20 random SNPs
  
    model1 = lm(Score ~SNP1+SNP2+SNP3+SNP4+SNP5+SNP6+SNP7+SNP8+SNP9+SNP10+SNP11+SNP12+SNP13+SNP14+SNP15+SNP16+SNP17+SNP18+SNP19+SNP20 ,data=merged_filt)
    PredictedValues = predict(model1)
    

    Actualmean = mean(na.omit(merged_filt$Score))
    Predictedmean= mean(PredictedValues)
  
    top = 0
    bottom = 0
  
    for (i in 1:nrow(merged_filt)){
    
      add_top = (PredictedValues[i]-Predictedmean)^2
      add_bottom = (merged_filt$Score[i]-Actualmean)^2
    
      top = top + add_top
      bottom = bottom + add_bottom
    }
  
    r2 = top /bottom
    PercentageVarianceExplained_20RandomMarkers = round(r2*100,1)
    PercentageVarianceExplained_20RandomMarkersSUM = PercentageVarianceExplained_20RandomMarkersSUM+PercentageVarianceExplained_20RandomMarkers
    
    if (PercentageVarianceExplained_20RandomMarkers>=PercentageVarianceExplained_AllMarkers){
      FDR_COUNT_20SelectionMarkers = FDR_COUNT_20SelectionMarkers+1
    }
  PVE_20_randomMarkers = PercentageVarianceExplained_20RandomMarkersSUM/1000
  Results_df[38] =PVE_20_randomMarkers
  FDR_20SelectionMarkers = FDR_COUNT_20SelectionMarkers/1000
  Results_df[39] =FDR_20SelectionMarkers
  
  if (exists('PercentageVarianceExplained_20TopGWASMarkers')){
    if (PercentageVarianceExplained_20RandomMarkers>=PercentageVarianceExplained_20TopGWASMarkers){
      FDR_COUNT_20GWASMarkers = FDR_COUNT_20GWASMarkers+1
    }
    FDR_20GWASMarkers = FDR_COUNT_20GWASMarkers/1000
    Results_df[42] =FDR_20GWASMarkers
  }
  }
  
  # Variance explained by EachLDBlock
  
  relevantSNPS = LDBlock1[c(1,4:7)] 
  #Then see how much all selection markers can explain
  Genotype_relevant <- Genotypes1[,c(1,which(colnames(Genotypes1)%in% relevantSNPS))]
  colnames(Genotype_relevant)[1]="Name"
  merged = merge(Phenotype_relevant,Genotype_relevant,by="Name",all = T)
  merged_filt = na.omit(merged)
  colnames(merged_filt)[1]=c("Line")
  merged_filt$Line=as.factor(merged_filt$Line)
  merged_filt$Score=as.numeric(as.character((merged_filt$Score)))
  
  colnames(merged_filt)[3:ncol(merged_filt)]=paste("SNP",seq(1:20),sep="")
  # Check variance explained by SNP
  
  model = lm(Score ~SNP1+SNP2+SNP3+SNP4+SNP5 ,data=merged_filt)
  
  PredictedValues = predict(model)
  
  Actualmean = mean(na.omit(merged_filt$Score))
  Predictedmean= mean(PredictedValues)
  
  top = 0
  bottom = 0
  
  for (i in 1:nrow(merged_filt)){
    
    add_top = (PredictedValues[i]-Predictedmean)^2
    add_bottom = (merged_filt$Score[i]-Actualmean)^2
    
    top = top + add_top
    bottom = bottom + add_bottom
  }
  
  r2 = top /bottom
  PercentageVarianceExplained_LDBlock1Markers = round(r2*100,1)
  
  Results_df[43] =PercentageVarianceExplained_LDBlock1Markers
  
  
  relevantSNPS = LDBlock2[c(1,4:5)] 
  #Then see how much all selection markers can explain
  Genotype_relevant <- Genotypes1[,c(1,which(colnames(Genotypes1)%in% relevantSNPS))]
  colnames(Genotype_relevant)[1]="Name"
  merged = merge(Phenotype_relevant,Genotype_relevant,by="Name",all = T)
  merged_filt = na.omit(merged)
  colnames(merged_filt)[1]=c("Line")
  merged_filt$Line=as.factor(merged_filt$Line)
  merged_filt$Score=as.numeric(as.character((merged_filt$Score)))
  
  colnames(merged_filt)[3:ncol(merged_filt)]=paste("SNP",seq(1:20),sep="")
  # Check variance explained by SNP
  
  model = lm(Score ~SNP1+SNP2+SNP3 ,data=merged_filt)
  
  PredictedValues = predict(model)
  
  Actualmean = mean(na.omit(merged_filt$Score))
  Predictedmean= mean(PredictedValues)
  
  top = 0
  bottom = 0
  
  for (i in 1:nrow(merged_filt)){
    
    add_top = (PredictedValues[i]-Predictedmean)^2
    add_bottom = (merged_filt$Score[i]-Actualmean)^2
    
    top = top + add_top
    bottom = bottom + add_bottom
  }
  
  r2 = top /bottom
  PercentageVarianceExplained_LDBlock2Markers = round(r2*100,1)
  
  Results_df[44] =PercentageVarianceExplained_LDBlock2Markers
  
  
  
  relevantSNPS = LDBlock3
  #Then see how much all selection markers can explain
  Genotype_relevant <- Genotypes1[,c(1,which(colnames(Genotypes1)%in% relevantSNPS))]
  colnames(Genotype_relevant)[1]="Name"
  merged = merge(Phenotype_relevant,Genotype_relevant,by="Name",all = T)
  merged_filt = na.omit(merged)
  colnames(merged_filt)[1]=c("Line")
  merged_filt$Line=as.factor(merged_filt$Line)
  merged_filt$Score=as.numeric(as.character((merged_filt$Score)))
  
  colnames(merged_filt)[3:ncol(merged_filt)]=paste("SNP",seq(1:20),sep="")
  # Check variance explained by SNP
  
  model = lm(Score ~SNP1+SNP2 ,data=merged_filt)
  
  PredictedValues = predict(model)
  
  Actualmean = mean(na.omit(merged_filt$Score))
  Predictedmean= mean(PredictedValues)
  
  top = 0
  bottom = 0
  
  for (i in 1:nrow(merged_filt)){
    
    add_top = (PredictedValues[i]-Predictedmean)^2
    add_bottom = (merged_filt$Score[i]-Actualmean)^2
    
    top = top + add_top
    bottom = bottom + add_bottom
  }
  
  r2 = top /bottom
  PercentageVarianceExplained_LDBlock3Markers = round(r2*100,1)
  
  Results_df[45] =PercentageVarianceExplained_LDBlock3Markers
  
  
  
  # 5 random markers that are not among LDBlock1 markers
    for (h in seq(1:1000)){
    
    SNPs_not_under_Selection <- colnames(Genotypes1)[-which(colnames(Genotypes1)%in% LDBlock1)]
    SNPs_not_under_Selection = SNPs_not_under_Selection[2:length(SNPs_not_under_Selection)]
    random_marker =sample(SNPs_not_under_Selection,5)
    Genotype_relevant <- Genotypes1[,c(1,which(colnames(Genotypes1)%in% random_marker))]
    colnames(Genotype_relevant)[1]="Name"
    merged = merge(Phenotype_relevant,Genotype_relevant,by="Name",all = T)
    merged_filt = na.omit(merged)
    colnames(merged_filt)[1]=c("Line")
    merged_filt$Line=as.factor(merged_filt$Line)
    merged_filt$Score=as.numeric(as.character((merged_filt$Score)))
    
    colnames(merged_filt)[3:ncol(merged_filt)]=paste("SNP",seq(1:20),sep="")
    # Check variance explained by all 20 random SNPs
    
    model1 = lm(Score ~SNP1+SNP2+SNP3+SNP4+SNP5 ,data=merged_filt)
    PredictedValues = predict(model1)
    
    
    Actualmean = mean(na.omit(merged_filt$Score))
    Predictedmean= mean(PredictedValues)
    
    top = 0
    bottom = 0
    
    for (i in 1:nrow(merged_filt)){
      
      add_top = (PredictedValues[i]-Predictedmean)^2
      add_bottom = (merged_filt$Score[i]-Actualmean)^2
      
      top = top + add_top
      bottom = bottom + add_bottom
    }
    
    r2 = top /bottom
    PercentageVarianceExplained_5RandomMarkers = round(r2*100,1)
    PercentageVarianceExplained_5RandomMarkersSUM = PercentageVarianceExplained_5RandomMarkersSUM+PercentageVarianceExplained_5RandomMarkers
    
    if (PercentageVarianceExplained_5RandomMarkers>=PercentageVarianceExplained_LDBlock1Markers){
      FDR_COUNT_5SelectionMarkers = FDR_COUNT_5SelectionMarkers+1
    }
    }
    PVE_5_randomMarkers = PercentageVarianceExplained_5RandomMarkersSUM/1000
    Results_df[46] =PVE_5_randomMarkers
    FDR_LDBlock1 = FDR_COUNT_5SelectionMarkers/1000
    Results_df[47] =FDR_LDBlock1
  
  
  # 3 random markers that are not among LDBlock1 markers
  for (h in seq(1:1000)){
    
    SNPs_not_under_Selection <- colnames(Genotypes1)[-which(colnames(Genotypes1)%in% LDBlock2)]
    SNPs_not_under_Selection = SNPs_not_under_Selection[2:length(SNPs_not_under_Selection)]
    random_marker =sample(SNPs_not_under_Selection,3)
    Genotype_relevant <- Genotypes1[,c(1,which(colnames(Genotypes1)%in% random_marker))]
    colnames(Genotype_relevant)[1]="Name"
    merged = merge(Phenotype_relevant,Genotype_relevant,by="Name",all = T)
    merged_filt = na.omit(merged)
    colnames(merged_filt)[1]=c("Line")
    merged_filt$Line=as.factor(merged_filt$Line)
    merged_filt$Score=as.numeric(as.character((merged_filt$Score)))
    
    colnames(merged_filt)[3:ncol(merged_filt)]=paste("SNP",seq(1:20),sep="")
    # Check variance explained by all 20 random SNPs
    
    model1 = lm(Score ~SNP1+SNP2+SNP3 ,data=merged_filt)
    PredictedValues = predict(model1)
    
    
    Actualmean = mean(na.omit(merged_filt$Score))
    Predictedmean= mean(PredictedValues)
    
    top = 0
    bottom = 0
    
    for (i in 1:nrow(merged_filt)){
      
      add_top = (PredictedValues[i]-Predictedmean)^2
      add_bottom = (merged_filt$Score[i]-Actualmean)^2
      
      top = top + add_top
      bottom = bottom + add_bottom
    }
    
    r2 = top /bottom
    PercentageVarianceExplained_3RandomMarkers = round(r2*100,1)
    PercentageVarianceExplained_3RandomMarkersSUM = PercentageVarianceExplained_3RandomMarkersSUM+PercentageVarianceExplained_3RandomMarkers
    
    if (PercentageVarianceExplained_3RandomMarkers>=PercentageVarianceExplained_LDBlock2Markers){
      FDR_COUNT_3SelectionMarkers = FDR_COUNT_3SelectionMarkers+1
    }
  }
    
    PVE_3_randomMarkers = PercentageVarianceExplained_3RandomMarkersSUM/1000
    Results_df[48] =PVE_3_randomMarkers
    FDR_LDBlock2 = FDR_COUNT_3SelectionMarkers/1000
    Results_df[49] =FDR_LDBlock2
  
  # 2 random markers that are not among LDBlock3 markers
  for (h in seq(1:1000)){
    
    SNPs_not_under_Selection <- colnames(Genotypes1)[-which(colnames(Genotypes1)%in% LDBlock3)]
    SNPs_not_under_Selection = SNPs_not_under_Selection[2:length(SNPs_not_under_Selection)]
    random_marker =sample(SNPs_not_under_Selection,2)
    Genotype_relevant <- Genotypes1[,c(1,which(colnames(Genotypes1)%in% random_marker))]
    colnames(Genotype_relevant)[1]="Name"
    merged = merge(Phenotype_relevant,Genotype_relevant,by="Name",all = T)
    merged_filt = na.omit(merged)
    colnames(merged_filt)[1]=c("Line")
    merged_filt$Line=as.factor(merged_filt$Line)
    merged_filt$Score=as.numeric(as.character((merged_filt$Score)))
    
    colnames(merged_filt)[3:ncol(merged_filt)]=paste("SNP",seq(1:20),sep="")
    # Check variance explained by all 20 random SNPs
    
    model1 = lm(Score ~SNP1+SNP2 ,data=merged_filt)
    PredictedValues = predict(model1)
    
    
    Actualmean = mean(na.omit(merged_filt$Score))
    Predictedmean= mean(PredictedValues)
    
    top = 0
    bottom = 0
    
    for (i in 1:nrow(merged_filt)){
      
      add_top = (PredictedValues[i]-Predictedmean)^2
      add_bottom = (merged_filt$Score[i]-Actualmean)^2
      
      top = top + add_top
      bottom = bottom + add_bottom
    }
    
    r2 = top /bottom
    PercentageVarianceExplained_2RandomMarkers = round(r2*100,1)
    PercentageVarianceExplained_2RandomMarkersSUM = PercentageVarianceExplained_2RandomMarkersSUM+PercentageVarianceExplained_2RandomMarkers
    
    if (PercentageVarianceExplained_2RandomMarkers>=PercentageVarianceExplained_LDBlock3Markers){
      FDR_COUNT_2SelectionMarkers = FDR_COUNT_2SelectionMarkers+1
    }
  }
    PVE_2_randomMarkers = PercentageVarianceExplained_2RandomMarkersSUM/1000
    Results_df[50] =PVE_2_randomMarkers
    FDR_LDBlock3 = FDR_COUNT_2SelectionMarkers/1000
    Results_df[51] =FDR_LDBlock3
  
  

  return(list(Results_df,Random_marker_df))
}



# Run function for all traits
for (trait in All_traits){
  print(trait)
  
  Temp_df1 <-PVE_function(trait,GWAS=1)[[1]]
  Temp_df2 <-PVE_function(trait,GWAS=1)[[2]]
  
  
  if (trait ==All_traits[1]){
    FullDF1 <- Temp_df1
    FullDF2 <- Temp_df2
    
  }else {
    FullDF1<-rbind(FullDF1, Temp_df1)
    FullDF2<-rbind(FullDF2, Temp_df2)
  }
  }


# Prepare for heatmap
rownames(FullDF1) = FullDF1[,1]
colnames(FullDF1)[(ncol(FullDF1)-14)]="All"
colnames(FullDF1)[(ncol(FullDF1)-13)]="20RandomMarkers"
colnames(FullDF1)[(ncol(FullDF1)-12)]="FDR_20SelectionMarkerss"
colnames(FullDF1)[(ncol(FullDF1)-11)]="Top20GWAS"
colnames(FullDF1)[(ncol(FullDF1)-10)]="1RandomMarker"
colnames(FullDF1)[(ncol(FullDF1)-9)]="FDR_Top20GWAS"
colnames(FullDF1)[(ncol(FullDF1)-8)]="LDBlock1"
colnames(FullDF1)[(ncol(FullDF1)-7)]="LDBlock2"
colnames(FullDF1)[(ncol(FullDF1)-6)]="LDBlock3"
colnames(FullDF1)[(ncol(FullDF1)-5)]="5RandomMarkers"
colnames(FullDF1)[(ncol(FullDF1)-4)]="FDR_LDBlock1"
colnames(FullDF1)[(ncol(FullDF1)-3)]="3RandomMarkers"
colnames(FullDF1)[(ncol(FullDF1)-2)]="FDR_LDBlock2"
colnames(FullDF1)[(ncol(FullDF1)-1)]="2RandomMarkers"
colnames(FullDF1)[(ncol(FullDF1))]="FDR_LDBlock3"


FullDF1_ <- FullDF1[,colSums(is.na(FullDF1))<nrow(FullDF1)]
FullDF1 = cbind(FullDF1_,FullDF1$Top20GWAS)


FullDF1_limited=FullDF1[,1:(ncol(FullDF1)-16)]
data <- data.matrix(FullDF1_limited)
col=colorRampPalette(rev(brewer.pal(n = 11, name =  "Spectral")))(1000)
p1 =pheatmap(data,color = col)
p1 # Full clustering present


# add suppopulation differentation
my_sample_col <- data.frame(SP = c("SP1","SP1","SP1","SP3","SP3","SP3","SP2_vs_SP3","SP3","SP3","SP3","SP3","SP3","SP3","SP2_vs_SP3","SP3","SP3","SP1","SP3","SP2","SP1","SP3","SP3","SP3","SP1","SP3"))
rownames(my_sample_col) <- colnames(FullDF1)[1:25]
my_colour = list(
  SP = c(SP1 = "#440154", SP2 = "#fde725",SP3="#21918c",SP2_vs_SP3="grey")
)


p1_alt1 =pheatmap::pheatmap(data, color=col,annotation_col = my_sample_col,annotation_colors = my_colour)
p1_alt1


RowAnnotation = data.frame(All= FullDF1$All ,LDBlock1=FullDF1$LDBlock1,LDBlock2 = FullDF1$LDBlock2,LDBlock3 = FullDF1$LDBlock3)
rownames(RowAnnotation) <- rownames(FullDF1)

my_colour = list(
  SP = c(SP1 = "#440154", SP2 = "#fde725",SP3="#21918c",SP2_vs_SP3="grey"),
  All = brewer.pal(n = 9, name = "YlOrBr"),
  LDBlock1 = brewer.pal(n = 9, name = "YlOrBr"),
  LDBlock2 = brewer.pal(n = 9, name = "YlOrBr"),
  LDBlock3 = brewer.pal(n = 9, name = "YlOrBr")
  )


p2 = pheatmap::pheatmap(data,color = col,cluster_cols=F)
p2 # Traits clustered, but markers in chromosomal order

data1=data[order(rownames(data)),]
p3 = pheatmap::pheatmap(data1,color = col,cluster_cols=F,cluster_rows=F,legend=T)
p3 # Traits in alphabetic order, markers in chromosomal order


# Save heatmap
plotname1 = paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/PVE_Heatmap_FullClustering_20markers_LDBlocks",Sys.Date(),".pdf",sep="")
plotname2 = paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/PVE_Heatmap_TraitClustering_20markers_LDBlocks",Sys.Date(),".pdf",sep="")
plotname3 = paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/PVE_Heatmap_NoClustering_20markers_LDBlocks",Sys.Date(),".pdf",sep="")
plotname4 = paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/PVE_Heatmap_FullClustering_Annotated_20markers_LDBlocks",Sys.Date(),".pdf",sep="")
plotname5 = paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/PVE_Heatmap_FullClustering_AnnotatedChrAsAnnotationGroup_20markers_LDBlocks",Sys.Date(),".pdf",sep="")
plotname6 = paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/PVE_Heatmap_FullClustering_AnnotatedChrAsAnnotationGroupWhite_20markers_LDBlocks",Sys.Date(),".pdf",sep="")
plotname7 = paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/PVE_Heatmap_FullClustering_AnnotatedChrAsAnnotationGroupWhite_randommarkers_20markers_LDBlocks",Sys.Date(),".pdf",sep="")
plotname8 = paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/PVE_Heatmap_FullClustering_AnnotatedChr_Significance1_20markers_LDBlocks",Sys.Date(),".pdf",sep="")
plotname9 = paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/PVE_Heatmap_FullClustering_AnnotatedChr_Significance2_20markers_LDBlocks",Sys.Date(),".pdf",sep="")
plotname10 = paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/PVE_Heatmap_FullClustering_AnnotatedCh_Significance3_20markers_LDBlocks",Sys.Date(),".pdf",sep="")
plotname11 = paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/PVE_Heatmap_FullClustering_AnnotatedCh_Significance3GWAS_20markers_LDBlocks",Sys.Date(),".pdf",sep="")
plotname12 = paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/PVE_Heatmap_FullClustering_AnnotatedCh_Significance3GWAS_20markers_LDBlocksCheating",Sys.Date(),".pdf",sep="")


save_pheatmap_pdf <- function(x, filename, width=10, height=10) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(p1, plotname1)
save_pheatmap_pdf(p2, plotname2)
save_pheatmap_pdf(p3, plotname3)

library(ComplexHeatmap)
pdf(plotname4, width=10, height=10)
p1_alt2 =ComplexHeatmap::pheatmap(data, color=col,annotation_col = my_sample_col,annotation_colors = my_colour,annotation_row = RowAnnotation,annotation_legend=TRUE)
p1_alt2
dev.off()


# Put annotation with chromosome colour as well
my_sample_col$Chr 
PositionsForMarkers <-SelMarkers[,c(1,2)]
colnames(PositionsForMarkers)=c("SNP","Chr")
my_sample_col$SNP=row.names(my_sample_col)
merged = merge(my_sample_col,PositionsForMarkers,by="SNP",sort=F)
rownames(merged) = merged[,1]
merged = merged[,-1]


my_colour = list(
  SP = c(SP1 = "#440154", SP2 = "#fde725",SP3="#21918c",SP2_vs_SP3="grey"),
  All = brewer.pal(n = 9, name = "YlOrBr"),
  LDBlock1 = brewer.pal(n = 9, name = "YlOrBr"),
  LDBlock2 = brewer.pal(n = 9, name = "YlOrBr"),
  LDBlock3 = brewer.pal(n = 9, name = "YlOrBr"),
  Chr = c(Chr1S="#b2182b",Chr1L="#ef8a62",Chr2="#fddbc7",Chr3="#f7f7f7",Chr4="#d1e5f0",Chr5="#67a9cf",Chr6="#2166ac"))

merged$Chr=as.factor(merged$Chr)
levels(merged$Chr)=c("Chr1S","Chr1L","Chr2","Chr3","Chr4","Chr5","Chr6")
pdf(plotname5, width=10, height=10)
p1_alt3 =ComplexHeatmap::pheatmap(data, color=col,annotation_col = merged,annotation_colors = my_colour,annotation_row = RowAnnotation,annotation_legend=TRUE)
p1_alt3
dev.off()

my_colour = list(
  SP = c(SP1 = "#440154", SP2 = "#fde725",SP3="#21918c",SP2_vs_SP3="grey"),
  All = brewer.pal(n = 9, name = "YlOrBr"),
  LDBlock1 = brewer.pal(n = 9, name = "YlOrBr"),
  LDBlock2 = brewer.pal(n = 9, name = "YlOrBr"),
  LDBlock3 = brewer.pal(n = 9, name = "YlOrBr"),
  Chr = c(Chr1S="white",Chr1L="white",Chr2="white",Chr3="white",Chr4="white",Chr5="white",Chr6="white"))

pdf(plotname6, width=10, height=10)
p1_alt3 =ComplexHeatmap::pheatmap(data, color=col,annotation_col = merged,annotation_colors = my_colour,annotation_row = RowAnnotation,annotation_legend=TRUE)
p1_alt3
dev.off()



# now add and asteriks to all cells that are significant at FDR  <0.005. That is less than 5/1000 times something random gives better results
mat2= data
mat2copy = data.matrix(mat2)

for (g in 1:nrow(data)){
  RandomMarkerPerc = FullDF2[g,2:ncol(FullDF2)]
  RandomMarkerPerc_sorted = RandomMarkerPerc[order(-RandomMarkerPerc)]
  FDR_cutoff_005 = as.numeric(RandomMarkerPerc_sorted[5])
  FDR_cutoff_01 = as.numeric(RandomMarkerPerc_sorted[10])
  FDR_cutoff_05 = as.numeric(RandomMarkerPerc_sorted[50])
  
  idx=which(mat2[g,1:20]>FDR_cutoff_005)
  if (length(idx)!=0){
    mat2copy[g,idx]="***"
  } 
  
  idx2=which(mat2[g,1:20]>FDR_cutoff_01  & mat2[g,1:20]<FDR_cutoff_005)
  if (length(idx2)!=0){
    mat2copy[g,idx2]="**"
  } 
  idx3=which(mat2[g,1:20]>FDR_cutoff_05 & mat2[g,1:20]<FDR_cutoff_01)
  if (length(idx3)!=0){
    mat2copy[g,idx3]="*"
  } 
  
  idx4=which(mat2[g,1:20]<FDR_cutoff_05)
  if (length(idx4)!=0){
    mat2copy[g,idx4]="NS"
  } 
}


pdf(plotname8, width=10, height=10)
p1_alt4 =ComplexHeatmap::pheatmap(data,color=col,annotation_col=merged,annotation_colors =my_colour,annotation_row=RowAnnotation,annotation_legend=TRUE, cell_fun = function(j, i, x, y, w, h, fill) {
  if(mat2copy[i, j] == "***") {
    grid.text("***", x, y)
  } 
  if(mat2copy[i, j] == "*") {
    grid.text("*", x, y)
  } 
  
  else if(mat2copy[i, j] =="**") {
    grid.text("**", x, y)
  }
})

p1_alt4
dev.off()



pdf(plotname9, width=10, height=10)
p1_alt5 =ComplexHeatmap::pheatmap(data,color=col,annotation_col=merged,annotation_colors =my_colour,annotation_row=RowAnnotation,annotation_legend=TRUE, cell_fun = function(j, i, x, y, w, h, fill) {
  if(mat2copy[i, j] == "***") {
    grid.text("***", x, y)
  } 
  
  else if(mat2copy[i, j] =="**") {
    grid.text("**", x, y)
  }
})

p1_alt5
dev.off()

pdf(plotname10, width=10, height=10)
p1_alt6 =ComplexHeatmap::pheatmap(data,color=col,annotation_col=merged,annotation_colors =my_colour,annotation_row=RowAnnotation,annotation_legend=TRUE, cell_fun = function(j, i, x, y, w, h, fill) {
  if(mat2copy[i, j] == "***") {
    grid.text("***", x, y)
  } 
})

p1_alt6
dev.off()



# And with GWAS markers

my_colour = list(
  SP = c(SP1 = "#440154", SP2 = "#fde725",SP3="#21918c",SP2_vs_SP3="grey"),
  All = brewer.pal(n = 9, name = "YlOrBr"),
  LDBlock1 = brewer.pal(n = 9, name = "YlOrBr"),
  LDBlock2 = brewer.pal(n = 9, name = "YlOrBr"),
  LDBlock3 = brewer.pal(n = 9, name = "YlOrBr"),
  Top20GWAS = brewer.pal(n = 9, name = "YlOrBr"),
  Chr = c(Chr1S="white",Chr1L="white",Chr2="white",Chr3="white",Chr4="white",Chr5="white",Chr6="white"))

RowAnnotation = data.frame(All=FullDF1$All,LDBlock1 = FullDF1$LDBlock1,LDBlock2 = FullDF1$LDBlock2,LDBlock3 = FullDF1$LDBlock3,Top20GWAS = FullDF1$Top20GWAS)
rownames(RowAnnotation) <- rownames(FullDF1)
RowAnnotation$Top20GWAS[is.na(RowAnnotation$Top20GWAS)]=0

pdf(plotname11, width=10, height=10)
p1_alt7 =ComplexHeatmap::pheatmap(data,color=col,annotation_col=merged,annotation_colors =my_colour,annotation_row=RowAnnotation,annotation_legend=TRUE, cell_fun = function(j, i, x, y, w, h, fill) {
  if(mat2copy[i, j] == "***") {
    grid.text("***", x, y)
  } 
  
  else if(mat2copy[i, j] =="**") {
    grid.text("**", x, y)
  }
})

p1_alt7
dev.off()


# a bit of cheating
cheatingvector1=(rep(max(RowAnnotation$Top20GWAS),5))
cheatingvector2=(rep(0,5))
RowAnnotation_new = rbind(RowAnnotation,cheatingvector1,cheatingvector2)
data1 =rep(0,ncol(data))
datacopy = rbind(data,data1,data1)
mat3copy = rbind(mat2copy,rep("NS",ncol(mat2copy)),rep("NS",ncol(mat2copy)))

pdf(plotname12, width=10, height=10)
p1_alt8 =ComplexHeatmap::pheatmap(datacopy,color=col,annotation_col=merged,annotation_colors =my_colour,annotation_row=RowAnnotation_new,annotation_legend=TRUE, cell_fun = function(j, i, x, y, w, h, fill) {
  if(mat3copy[i, j] == "***") {
    grid.text("***", x, y)
  } 
  
  else if(mat3copy[i, j] =="**") {
    grid.text("**", x, y)
  }
})

p1_alt8
dev.off()

# Save file with PVE of each trait pr. SNP
filename= paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/PVE_All_traits_SelectionMarkers_20Markers_LDBlocks",Sys.Date(),".txt",sep="")
write.table(FullDF1,filename,sep="\t",quote=F,row.names = T,col.names = T)

# Save file with PVE of each trait pr. random SNP
filename2= paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/PVE_All_traits_OneRandomMarker_20Markers_LDBlocks",Sys.Date(),".txt",sep="")
write.table(FullDF2,filename2,sep="\t",quote=F,row.names = T,col.names = T)













# Check if single marker explains more than expected by random

SignificantPVEMarkers<- data.frame(matrix(ncol = 3, nrow = 1))
colnames(SignificantPVEMarkers) = c("trait","Marker","FDR")

for (q in seq(1:nrow(FullDF1))){
  
  trait = rownames(FullDF1)[q]
  for (h in seq(1:20)){
    givenmarker = colnames(FullDF1)[h]
    givenmarker_PVE = FullDF1[q,h]
    FDR_For_Marker = length(which(FullDF2[q,2:ncol(FullDF2)]>=givenmarker_PVE))/1000
    
    if (FDR_For_Marker<0.05){
      print(paste(givenmarker,"is significant for trait:",trait,"FDR=",FDR_For_Marker,sep=" "))
      SignificantPVEMarkers_temp = c(trait,givenmarker,FDR_For_Marker)
      SignificantPVEMarkers = rbind(SignificantPVEMarkers,SignificantPVEMarkers_temp)
      
    }
  }
}
filename3= paste("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220425_PVE_SelectionMarkers_Heatmap/SignificantMarkers_20Markers_LDblocks",Sys.Date(),".txt",sep="")
write.table(SignificantPVEMarkers,filename3,sep="\t",quote=F,row.names = T,col.names = T)



