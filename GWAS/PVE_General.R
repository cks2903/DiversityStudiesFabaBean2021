#Phenotypes = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/GWAS/Pheno.csv",sep=",",head=T)
Phenotypes = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/_Phenotypes_qualityfiltered_minimum_2022-03-25_.csv",sep=",",head=T)
Phenotypes$PDIDTRID = paste(Phenotypes$PDID,Phenotypes$TRID,sep=":")

Genotypes1 = read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220324_GWAS_MAGIC/GWAS/Genotypes_FarmCPU.txt",sep="\t",head=T, check.names = F)
colnames(Genotypes1)


PVE_function <- function(trait,TRID=NULL, SNP){
  
  if (length(TRID)==0){
    Phenotype_relevant <- Phenotypes[which(Phenotypes$PDID==trait),]
    # average over lines
    #Phenotype_relevant <- aggregate(Phenotype_relevant[, 3], list(Phenotype_relevant$Name), mean)

  } else{
    Phenotype_relevant <- Phenotypes[which(Phenotypes$PDID==trait & Phenotypes$TRID==TRID),]
  }
  Phenotype_relevant =Phenotype_relevant[,c(8,3)]
  colnames(Phenotype_relevant)=c("Name","Score")
  Genotype_relevant <- Genotypes1[,c(1,which(colnames(Genotypes1)==SNP))]
  colnames(Genotype_relevant)[1]="Name"
    
    merged = merge(Phenotype_relevant,Genotype_relevant,by="Name",all = T)
    merged_filt = na.omit(merged)
    colnames(merged_filt)=c("Line","Score","Marker")
    merged_filt$Marker=as.factor(merged_filt$Marker)
    merged_filt$Line=as.factor(merged_filt$Line)
    merged_filt$Score=as.numeric(as.character((merged_filt$Score)))
    
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
    print(r2*100)
    
    
  }
  

# Examples of how to use the function
PVE_function(trait="29", TRID=NULL,SNP="AX-416788629")
PVE_function(trait="29", TRID=NULL,SNP="AX-181482207")
PVE_function(trait="29", TRID=NULL,SNP="AX-416797900")


PVE_function(trait="27", TRID=NULL,SNP="AX-416816373")
PVE_function(trait="27", TRID=NULL,SNP="AX-181482848")
PVE_function(trait="27", TRID=NULL,SNP="AX-416763140")
PVE_function(trait="27", TRID=NULL,SNP="AX-416742669")
PVE_function(trait="27", TRID=NULL,SNP="AX-181153730")


PVE_function(trait="27", TRID=31,SNP="AX-416729972")
PVE_function(trait="27", TRID=31,SNP="AX-181185226")
PVE_function(trait="27", TRID=31,SNP="AX-181489476")
PVE_function(trait="27", TRID=31,SNP="AX-416770212")
PVE_function(trait="27", TRID=31,SNP="AX-416758256")
PVE_function(trait="27", TRID=31,SNP="AX-416814584")
PVE_function(trait="27", TRID=31,SNP="AX-416825912")


PVE_function(trait="29", TRID=27,SNP="AX-416760427")
PVE_function(trait="29", TRID=27,SNP="AX-181185889")
PVE_function(trait="29", TRID=27,SNP="AX-416734975")
PVE_function(trait="29", TRID=27,SNP="AX-181177377")
PVE_function(trait="29", TRID=27,SNP="AX-416759210")


PVE_function(trait="28", TRID=c(28.2,28.1),SNP="AX-181200574")
PVE_function(trait="28", TRID=c(28.2,28.1),SNP="AX-416774097")
PVE_function(trait="28", TRID=c(28,1,28.2),SNP="AX-416758489")


PVE_function(trait="36", TRID=NULL,SNP="AX-416819177")
PVE_function(trait="36", TRID=NULL,SNP="AX-181159680")
PVE_function(trait="36", TRID=NULL,SNP="AX-416747244")
PVE_function(trait="36", TRID=NULL,SNP="AX-416820991")
PVE_function(trait="36", TRID=NULL,SNP="AX-416728514")
PVE_function(trait="36", TRID=NULL,SNP="AX-181460581")
PVE_function(trait="36", TRID=NULL,SNP="AX-181193226")

PVE_function(trait="36", TRID=31,SNP="AX-181201469")
PVE_function(trait="36", TRID=31,SNP="AX-416733502")
PVE_function(trait="36", TRID=31,SNP="AX-416747244")
PVE_function(trait="36", TRID=31,SNP="AX-181189191")
PVE_function(trait="36", TRID=31,SNP="AX-416738106")
PVE_function(trait="36", TRID=31,SNP="AX-181175084")

PVE_function(trait="36", TRID=29,SNP="AX-181485131")
PVE_function(trait="36", TRID=29,SNP="AX-416747244")
PVE_function(trait="36", TRID=29,SNP="AX-181159943")
PVE_function(trait="36", TRID=29,SNP="AX-416727937")
PVE_function(trait="36", TRID=29,SNP="AX-416746325")

PVE_function(trait="46", TRID=NULL,SNP="AX-181484321")
PVE_function(trait="46", TRID=NULL,SNP="AX-181463708")
PVE_function(trait="46", TRID=NULL,SNP="AX-416817171")
PVE_function(trait="46", TRID=NULL,SNP="AX-181189191")
PVE_function(trait="46", TRID=NULL,SNP="AX-181191536")
PVE_function(trait="46", TRID=NULL,SNP="AX-416735217")
PVE_function(trait="46", TRID=NULL,SNP="AX-181175083")

PVE_function(trait="46", TRID=28.1,SNP="AX-181484321")
PVE_function(trait="46", TRID=28.1,SNP="AX-416733812")
PVE_function(trait="46", TRID=28.1,SNP="AX-416817171")
PVE_function(trait="46", TRID=28.1,SNP="AX-181189191")

PVE_function(trait="46", TRID=28.2,SNP="AX-181488110")
PVE_function(trait="46", TRID=28.2,SNP="AX-416763828")
PVE_function(trait="46", TRID=28.2,SNP="AX-181189191")
PVE_function(trait="46", TRID=28.2,SNP="AX-181191536")
PVE_function(trait="46", TRID=28.2,SNP="AX-416760328")
PVE_function(trait="46", TRID=28.2,SNP="AX-416735217")

PVE_function(trait="53", TRID=NULL,SNP="AX-416807336")
PVE_function(trait="53", TRID=NULL,SNP="AX-416736377")
PVE_function(trait="53", TRID=NULL,SNP="AX-416753984")
PVE_function(trait="53", TRID=NULL,SNP="AX-416813947")
PVE_function(trait="53", TRID=NULL,SNP="AX-416796127")


PVE_function(trait="14", TRID=NULL,SNP="AX-416796020")
PVE_function(trait="14", TRID=NULL,SNP="AX-181446054")
PVE_function(trait="14", TRID=NULL,SNP="AX-181187562")
PVE_function(trait="14", TRID=NULL,SNP="AX-416787239")
PVE_function(trait="14", TRID=NULL,SNP="AX-416813291")
PVE_function(trait="14", TRID=NULL,SNP="AX-416777072")


PVE_function(trait="14", TRID=28,SNP="AX-416772711")
PVE_function(trait="14", TRID=28,SNP="AX-181490582")
PVE_function(trait="14", TRID=28,SNP="AX-416741468")
PVE_function(trait="14", TRID=28,SNP="AX-416800906")
PVE_function(trait="14", TRID=28,SNP="AX-416819455")

PVE_function(trait="14", TRID=31,SNP="AX-416796988")
PVE_function(trait="14", TRID=31,SNP="AX-416751143")
PVE_function(trait="14", TRID=31,SNP="AX-416806084")
PVE_function(trait="14", TRID=31,SNP="AX-181187562")
PVE_function(trait="14", TRID=31,SNP="AX-181488372")
PVE_function(trait="14", TRID=31,SNP="AX-416745649")
PVE_function(trait="14", TRID=31,SNP="AX-416791529")

PVE_function(trait="38", TRID=28,SNP="AX-416796020")
PVE_function(trait="38", TRID=28,SNP="AX-181445488")
PVE_function(trait="38", TRID=28,SNP="AX-181152545")


PVE_function(trait="40", TRID=NULL,SNP="AX-416775708")
PVE_function(trait="40", TRID=NULL,SNP="AX-181197227")
PVE_function(trait="40", TRID=NULL,SNP="AX-416757116")
PVE_function(trait="40", TRID=NULL,SNP="AX-416803185")

PVE_function(trait="40", TRID=28,SNP="AX-181464360")
PVE_function(trait="40", TRID=28,SNP="AX-181197227")
PVE_function(trait="40", TRID=28,SNP="AX-416791921")
PVE_function(trait="40", TRID=28,SNP="AX-416745463")
PVE_function(trait="40", TRID=28,SNP="AX-416795821")
PVE_function(trait="40", TRID=28,SNP="AX-416743059")

PVE_function(trait="17", TRID=NULL,SNP="AX-416722153")
PVE_function(trait="17", TRID=NULL,SNP="AX-416722950")
PVE_function(trait="17", TRID=NULL,SNP="AX-416824308")
PVE_function(trait="17", TRID=NULL,SNP="AX-181153171")

PVE_function(trait="17", TRID=28,SNP="AX-181487046")
PVE_function(trait="17", TRID=28,SNP="AX-416768891")
PVE_function(trait="17", TRID=28,SNP="AX-416722950")
PVE_function(trait="17", TRID=28,SNP="AX-181153171")
PVE_function(trait="17", TRID=28,SNP="AX-416766046")
PVE_function(trait="17", TRID=28,SNP="AX-416740003")
PVE_function(trait="17", TRID=28,SNP="AX-416732452")

PVE_function(trait="17", TRID=31,SNP="AX-416756553")
PVE_function(trait="17", TRID=31,SNP="AX-416768891")
PVE_function(trait="17", TRID=31,SNP="AX-416815482")
PVE_function(trait="17", TRID=31,SNP="AX-181153171")
PVE_function(trait="17", TRID=31,SNP="AX-416723283")
PVE_function(trait="17", TRID=31,SNP="AX-416787239")
PVE_function(trait="17", TRID=31,SNP="AX-181194859")

PVE_function(trait="17", TRID=29,SNP="AX-416787525")
PVE_function(trait="17", TRID=29,SNP="AX-181185917")
PVE_function(trait="17", TRID=29,SNP="AX-416723283")


PVE_function(trait="19", TRID=NULL,SNP="AX-416754977")
PVE_function(trait="19", TRID=NULL,SNP="AX-416726585")
PVE_function(trait="19", TRID=NULL,SNP="AX-181462618")
PVE_function(trait="19", TRID=NULL,SNP="AX-181483910")
PVE_function(trait="19", TRID=NULL,SNP="AX-416747467")
PVE_function(trait="19", TRID=NULL,SNP="AX-416811554")

PVE_function(trait="19", TRID=28,SNP="AX-181177649")
PVE_function(trait="19", TRID=28,SNP="AX-416810414")
PVE_function(trait="19", TRID=28,SNP="AX-416790910")
PVE_function(trait="19", TRID=28,SNP="AX-416759775")
PVE_function(trait="19", TRID=28,SNP="AX-416721892")

PVE_function(trait="19", TRID=27,SNP="AX-416763985")
PVE_function(trait="19", TRID=27,SNP="AX-181481959")
PVE_function(trait="19", TRID=27,SNP="AX-416732539")
PVE_function(trait="19", TRID=27,SNP="AX-416726585")
PVE_function(trait="19", TRID=27,SNP="AX-416815578")
PVE_function(trait="19", TRID=27,SNP="AX-416722630")
PVE_function(trait="19", TRID=27,SNP="AX-181486012")
PVE_function(trait="19", TRID=27,SNP="AX-181490096")
PVE_function(trait="19", TRID=27,SNP="AX-181168340")

PVE_function(trait="19", TRID=31,SNP="AX-181168340")
PVE_function(trait="19", TRID=31,SNP="AX-416726585")
PVE_function(trait="19", TRID=31,SNP="AX-416744502")
PVE_function(trait="19", TRID=31,SNP="AX-181464345")


PVE_function(trait="19", TRID=29,SNP="AX-181178634")
PVE_function(trait="19", TRID=29,SNP="AX-181481959")
PVE_function(trait="19", TRID=29,SNP="AX-416726585")
PVE_function(trait="19", TRID=29,SNP="AX-181150007")
PVE_function(trait="19", TRID=29,SNP="AX-181483910")
PVE_function(trait="19", TRID=29,SNP="AX-416760644")


PVE_function(trait="20", TRID=NULL,SNP="AX-181481959")
PVE_function(trait="20", TRID=NULL,SNP="AX-416722749")
PVE_function(trait="20", TRID=NULL,SNP="AX-181487107")
PVE_function(trait="20", TRID=NULL,SNP="AX-181193698")
PVE_function(trait="20", TRID=NULL,SNP="AX-181487700")
PVE_function(trait="20", TRID=NULL,SNP="AX-181483910")
PVE_function(trait="20", TRID=NULL,SNP="AX-181182248")


PVE_function(trait="20", TRID=28,SNP="AX-416821697")
PVE_function(trait="20", TRID=28,SNP="AX-416804927")
PVE_function(trait="20", TRID=28,SNP="AX-416722749")
PVE_function(trait="20", TRID=28,SNP="AX-181482425")
PVE_function(trait="20", TRID=28,SNP="AX-416737202")
PVE_function(trait="20", TRID=28,SNP="AX-181178810")
PVE_function(trait="20", TRID=28,SNP="AX-416721892")

PVE_function(trait="20", TRID=27,SNP="AX-416791615")
PVE_function(trait="20", TRID=27,SNP="AX-416737202")
PVE_function(trait="20", TRID=27,SNP="AX-181178810")
PVE_function(trait="20", TRID=27,SNP="AX-416744933")
PVE_function(trait="20", TRID=27,SNP="AX-181485141")


PVE_function(trait="20", TRID=31,SNP="AX-181180062")
PVE_function(trait="20", TRID=31,SNP="AX-416759775")
PVE_function(trait="20", TRID=31,SNP="AX-416730523")
PVE_function(trait="20", TRID=31,SNP="AX-416817327")
PVE_function(trait="20", TRID=31,SNP="AX-181487700")
PVE_function(trait="20", TRID=31,SNP="AX-416822997")
PVE_function(trait="20", TRID=31,SNP="AX-416749838")

PVE_function(trait="20", TRID=29,SNP="AX-181439302")
PVE_function(trait="20", TRID=29,SNP="AX-416733812")
PVE_function(trait="20", TRID=29,SNP="AX-416744502")

PVE_function(trait="22", TRID=NULL,SNP="AX-416754977")
PVE_function(trait="22", TRID=NULL,SNP="AX-181194033")
PVE_function(trait="22", TRID=NULL,SNP="AX-416780606")
PVE_function(trait="22", TRID=NULL,SNP="AX-181178807")
PVE_function(trait="22", TRID=NULL,SNP="AX-181487700")
PVE_function(trait="22", TRID=NULL,SNP="AX-416747267")
PVE_function(trait="22", TRID=NULL,SNP="AX-416789448")

PVE_function(trait="22", TRID=28,SNP="AX-416792689")
PVE_function(trait="22", TRID=28,SNP="AX-416751143")
PVE_function(trait="22", TRID=28,SNP="AX-416780606")
PVE_function(trait="22", TRID=28,SNP="AX-181446054")
PVE_function(trait="22", TRID=28,SNP="AX-181474476")
PVE_function(trait="22", TRID=28,SNP="AX-416776148")
PVE_function(trait="22", TRID=28,SNP="AX-181171286")
PVE_function(trait="22", TRID=28,SNP="AX-181487700")
PVE_function(trait="22", TRID=28,SNP="AX-181149099")

PVE_function(trait="22", TRID=31,SNP="AX-181482671")
PVE_function(trait="22", TRID=31,SNP="AX-416811006")

PVE_function(trait="22", TRID=31,SNP="AX-416810350")
PVE_function(trait="22", TRID=31,SNP="AX-416767693")
PVE_function(trait="22", TRID=31,SNP="AX-181446054")
PVE_function(trait="22", TRID=31,SNP="AX-416751954")
PVE_function(trait="22", TRID=31,SNP="AX-416755850")
PVE_function(trait="22", TRID=31,SNP="AX-181464345")

PVE_function(trait="21", TRID=NULL,SNP="AX-181481959")
PVE_function(trait="21", TRID=NULL,SNP="AX-181170911")
PVE_function(trait="21", TRID=NULL,SNP="AX-181487107")
PVE_function(trait="21", TRID=NULL,SNP="AX-416796690")
PVE_function(trait="21", TRID=NULL,SNP="AX-181193698")
PVE_function(trait="21", TRID=NULL,SNP="AX-181487700")
PVE_function(trait="21", TRID=NULL,SNP="AX-416814129")
PVE_function(trait="21", TRID=NULL,SNP="AX-181483910")
PVE_function(trait="21", TRID=NULL,SNP="AX-181182248")

PVE_function(trait="21", TRID=28,SNP="AX-181154521")
PVE_function(trait="21", TRID=28,SNP="AX-416804927")
PVE_function(trait="21", TRID=28,SNP="AX-416766902")
PVE_function(trait="21", TRID=28,SNP="AX-416785709")
PVE_function(trait="21", TRID=28,SNP="AX-181178810")
PVE_function(trait="21", TRID=28,SNP="AX-416722862")
PVE_function(trait="21", TRID=28,SNP="AX-416744931")
PVE_function(trait="21", TRID=28,SNP="AX-416753078")
PVE_function(trait="21", TRID=28,SNP="AX-181484352")

PVE_function(trait="21", TRID=27,SNP="AX-416737202")
PVE_function(trait="21", TRID=27,SNP="AX-181490096")
PVE_function(trait="21", TRID=27,SNP="AX-416722350")

PVE_function(trait="21", TRID=31,SNP="AX-181200574")
PVE_function(trait="21", TRID=31,SNP="AX-181173675")
PVE_function(trait="21", TRID=31,SNP="AX-416813841")
PVE_function(trait="21", TRID=31,SNP="AX-416726585")
PVE_function(trait="21", TRID=31,SNP="AX-416796690")
PVE_function(trait="21", TRID=31,SNP="AX-416755850")


PVE_function(trait="21", TRID=29,SNP="AX-416726585")
PVE_function(trait="21", TRID=29,SNP="AX-181439302")
PVE_function(trait="21", TRID=29,SNP="AX-416785530")
PVE_function(trait="21", TRID=29,SNP="AX-416744502")
PVE_function(trait="21", TRID=29,SNP="AX-181464345")
PVE_function(trait="21", TRID=29,SNP="AX-416729988")
PVE_function(trait="21", TRID=29,SNP="AX-181483910")
PVE_function(trait="21", TRID=29,SNP="AX-416762454")
PVE_function(trait="21", TRID=29,SNP="AX-416804743")


