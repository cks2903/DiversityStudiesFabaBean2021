library(varhandle)
library(tidyverse)
library(ggplot2)
library(viridis)
library(plyr)





####################################################
####################################################
######                K = 3                   ###### 
####################################################
####################################################

tbl=read.table("genotypesADMIXTURE_MAFfiltered.3.Q")


plot_data <- tbl %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V3) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
                assingment_prob = max(prob)) %>% 
  dplyr::arrange(likely_assignment, desc(assingment_prob)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(id = forcats::fct_inorder(factor(id)))

p3 = ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c")) +
  theme_classic()

p3
ggsave("Admixture_Vir_K=3.pdf", plot = p3, width = 40, height = 10, unit = 'cm')



# Plotting K=3 after origin instead
origin <- read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220224_ADMIXTURE/Origin_Broad.txt",head=T,sep="\t") #Load file with origin of different samplles
indviduals <- read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220224_ADMIXTURE/BedFileForStructureFiltered_20220223.fam",head=F,sep=" ")

tbl_ind = cbind(indviduals[,1],tbl)
colnames(tbl_ind) = c("Sample","Pop1","Pop2","Pop3")

Merged_origin <- join(tbl_ind,origin)
all(tbl_ind[,1]==Merged_origin[,1])

# Asia
Asia_Samples =Merged_origin[which(Merged_origin$OriginLessBroad=="Asia"),1]
tbl_asia = tbl[which(tbl_ind[,1] %in% Asia_Samples),]

plot_data_Asia <- tbl_asia %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V3) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
                assingment_prob = max(prob)) %>% 
  dplyr::arrange(likely_assignment, desc(assingment_prob)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(id = forcats::fct_inorder(factor(id)))

plot_data_Asia$order = rep("NA",nrow(plot_data_Asia))
count=1
for (i in seq(1,length(levels(plot_data_Asia$id)))){
  old_level = levels(plot_data_Asia$id)[i]
  idx=which(plot_data_Asia$id==old_level)
  plot_data_Asia$order[idx]=count
  count=count+1
}
plot_data_Asia$order=as.factor(as.integer(as.character(plot_data_Asia$order)))

pa = ggplot(plot_data_Asia, aes(order, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c")) +
  theme_classic()

pa
ggsave("Admixture_Asia_MiddleEastIncluded_K=3.pdf", plot = pa, width = 40, height = 10, unit = 'cm')


# North
North =Merged_origin[which(Merged_origin$OriginLessBroad=="North"),1]
tbl_North= tbl[which(tbl_ind[,1] %in% North),]



plot_data_North <- tbl_North %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V3) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
                assingment_prob = max(prob)) %>% 
  dplyr::arrange(likely_assignment, desc(assingment_prob)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(id = forcats::fct_inorder(factor(id)))


plot_data_North$order = rep("NA",nrow(plot_data_North))
count=1
for (i in seq(1,length(levels(plot_data_North$id)))){
  old_level = levels(plot_data_North$id)[i]
  idx=which(plot_data_North$id==old_level)
  plot_data_North$order[idx]=count
  count=count+1
}
plot_data_North$order=as.factor(as.integer(as.character(plot_data_North$order)))


pn = ggplot(plot_data_North, aes(order, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c")) +
  theme_classic()

pn
ggsave("Admixture_MiddleEastIncluded_North_K=3.pdf", plot = pn, width = 40, height = 10, unit = 'cm')



# South
South_Samples =Merged_origin[which(Merged_origin$OriginLessBroad=="South"),1]
tbl_South = tbl[which(tbl_ind[,1] %in% South_Samples),]



plot_data_South <- tbl_South %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V3) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
                assingment_prob = max(prob)) %>% 
  dplyr::arrange(likely_assignment, desc(assingment_prob)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(id = forcats::fct_inorder(factor(id)))



plot_data_South$order = rep("NA",nrow(plot_data_South))
count=1
for (i in seq(1,length(levels(plot_data_South$id)))){
  old_level = levels(plot_data_South$id)[i]
  idx=which(plot_data_South$id==old_level)
  plot_data_South$order[idx]=count
  count=count+1
}
plot_data_South$order=as.factor(as.integer(as.character(plot_data_South$order)))


ps = ggplot(plot_data_South, aes(order, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c")) +
  theme_classic()

ps
ggsave("Admixture_MiddleEastIncluded_South_K=3.pdf", plot = ps, width = 40, height = 10, unit = 'cm')



# MiddleEast
ME_Samples =Merged_origin[which(Merged_origin$OriginLessBroad=="Middle East"),1]
tbl_ME = tbl[which(tbl_ind[,1] %in% ME_Samples),]



plot_data_Me <- tbl_ME %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V3) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
                assingment_prob = max(prob)) %>% 
  dplyr::arrange(likely_assignment, desc(assingment_prob)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(id = forcats::fct_inorder(factor(id)))



plot_data_Me$order = rep("NA",nrow(plot_data_Me))
count=1
for (i in seq(1,length(levels(plot_data_Me$id)))){
  old_level = levels(plot_data_Me$id)[i]
  idx=which(plot_data_Me$id==old_level)
  plot_data_Me$order[idx]=count
  count=count+1
}
plot_data_Me$order=as.factor(as.integer(as.character(plot_data_Me$order)))


pme = ggplot(plot_data_Me, aes(order, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c")) +
  theme_classic()

pme
ggsave("Admixture_MiddleEastIncluded_MiddleEast_K=3.pdf", plot = pme, width = 40, height = 10, unit = 'cm')


# No location
NoLocation =Merged_origin[which(is.na(Merged_origin$OriginLessBroad=="South")),1]
tbl_NoLoc= tbl[which(tbl_ind[,1] %in% NoLocationIdx1),]

NoLocation_Samples = tbl[which(tbl_ind[,1] %in% NoLocation),]

plot_data_NoLoc <- NoLocation_Samples %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V3) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
                assingment_prob = max(prob)) %>% 
  dplyr::arrange(likely_assignment, desc(assingment_prob)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(id = forcats::fct_inorder(factor(id)))

plot_data_NoLoc$order = rep("NA",nrow(plot_data_NoLoc))
count=1
for (i in seq(1,length(levels(plot_data_NoLoc$id)))){
  old_level = levels(plot_data_NoLoc$id)[i]
  idx=which(plot_data_NoLoc$id==old_level)
  plot_data_NoLoc$order[idx]=count
  count=count+1
}
plot_data_NoLoc$order=as.factor(as.integer(as.character(plot_data_NoLoc$order)))


pnl = ggplot(plot_data_NoLoc, aes(order, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c")) +
  theme_classic()

pnl
ggsave("Admixture_MiddleEastIncluded_NoLocation_K=3.pdf", plot = pnl, width = 40, height = 10, unit = 'cm')





# bring together
Mixed_3 = rbind(plot_data_North,plot_data_South,plot_data_Asia,plot_data_Me,plot_data_NoLoc)
levels(Mixed_3$order)
Mixed_3$Location = rep("NA",nrow(Mixed_3))
Mixed_3$Location[1:420]="North"
Mixed_3$Location[421:(492+420)]="South"
Mixed_3$Location[913:(219+912)]="Asia"
Mixed_3$Location[1132:1218]="Middle East"


Mixed_3$order = unfactor(Mixed_3$order)
Mixed_3$order[421:(492+420)] = Mixed_3$order[421:(492+420)] +1000 # larger level numbers for South
Mixed_3$order[913:(219+912)] = Mixed_3$order[913:(219+912)] +100000 # even larger level numbers for Asia
Mixed_3$order[1132:1218] = Mixed_3$order[1132:1218] +10000 # even even larger level numbers for Middle east
Mixed_3$order = as.factor(Mixed_3$order)

Mixed_3$Location = as.factor(Mixed_3$Location)
Mixed_3$Location <- factor(Mixed_3$Location, levels = c("North", "South", "Middle East","Asia","NA"))

K3AllGeo = ggplot(Mixed_3, aes(order, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c")) +
  theme_classic()  +facet_grid(. ~ Location,scales = "free", space = "free_x")  +  theme(axis.text.x = element_text(angle = 90))


ggsave("Admixture_ByGeographicLocation_K=3.pdf", plot = K3AllGeo, width = 100, height = 20, unit = 'cm')











####################################################
####################################################
######                K = 4                   ###### 
####################################################
####################################################

tbl=read.table("genotypesADMIXTURE_MAFfiltered.4.Q")


plot_data <- tbl %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V4) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
                assingment_prob = max(prob)) %>% 
  dplyr::arrange(likely_assignment, desc(assingment_prob)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(id = forcats::fct_inorder(factor(id)))

p4 = ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c","#5ec962")) +
  theme_classic()

p4
ggsave("Admixture_Vir_K=4.pdf", plot = p4, width = 40, height = 10, unit = 'cm')



# Plotting K=4 after origin instead
origin <- read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220224_ADMIXTURE/Origin_Broad.txt",head=T,sep="\t")
indviduals <- read.table("/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220224_ADMIXTURE/BedFileForStructureFiltered_20220223.fam",head=F,sep=" ")

tbl_ind = cbind(indviduals[,1],tbl)
colnames(tbl_ind) = c("Sample","Pop1","Pop2","Pop3","Pop4")

Merged_origin <- join(tbl_ind,origin)
all(tbl_ind[,1]==Merged_origin[,1])

# Asia
Asia_Samples =Merged_origin[which(Merged_origin$OriginLessBroad=="Asia"),1]
tbl_asia = tbl[which(tbl_ind[,1] %in% Asia_Samples),]

plot_data_Asia <- tbl_asia %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V4) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
                assingment_prob = max(prob)) %>% 
  dplyr::arrange(likely_assignment, desc(assingment_prob)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(id = forcats::fct_inorder(factor(id)))

plot_data_Asia$order = rep("NA",nrow(plot_data_Asia))
count=1
for (i in seq(1,length(levels(plot_data_Asia$id)))){
  old_level = levels(plot_data_Asia$id)[i]
  idx=which(plot_data_Asia$id==old_level)
  plot_data_Asia$order[idx]=count
  count=count+1
}
plot_data_Asia$order=as.factor(as.integer(as.character(plot_data_Asia$order)))

pa = ggplot(plot_data_Asia, aes(order, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c","#5ec962")) +
  theme_classic()

pa
ggsave("Admixture_Asia_MiddleEastIncluded_K=4.pdf", plot = pa, width = 40, height = 10, unit = 'cm')


# North
North =Merged_origin[which(Merged_origin$OriginLessBroad=="North"),1]
tbl_North= tbl[which(tbl_ind[,1] %in% North),]



plot_data_North <- tbl_North %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V4) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
                assingment_prob = max(prob)) %>% 
  dplyr::arrange(likely_assignment, desc(assingment_prob)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(id = forcats::fct_inorder(factor(id)))


plot_data_North$order = rep("NA",nrow(plot_data_North))
count=1
for (i in seq(1,length(levels(plot_data_North$id)))){
  old_level = levels(plot_data_North$id)[i]
  idx=which(plot_data_North$id==old_level)
  plot_data_North$order[idx]=count
  count=count+1
}
plot_data_North$order=as.factor(as.integer(as.character(plot_data_North$order)))


pn = ggplot(plot_data_North, aes(order, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c","#5ec962")) +
  theme_classic()

pn
ggsave("Admixture_MiddleEastIncluded_North_K=4.pdf", plot = pn, width = 40, height = 10, unit = 'cm')



# South
South_Samples =Merged_origin[which(Merged_origin$OriginLessBroad=="South"),1]
tbl_South = tbl[which(tbl_ind[,1] %in% South_Samples),]



plot_data_South <- tbl_South %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V4) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
                assingment_prob = max(prob)) %>% 
  dplyr::arrange(likely_assignment, desc(assingment_prob)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(id = forcats::fct_inorder(factor(id)))



plot_data_South$order = rep("NA",nrow(plot_data_South))
count=1
for (i in seq(1,length(levels(plot_data_South$id)))){
  old_level = levels(plot_data_South$id)[i]
  idx=which(plot_data_South$id==old_level)
  plot_data_South$order[idx]=count
  count=count+1
}
plot_data_South$order=as.factor(as.integer(as.character(plot_data_South$order)))


ps = ggplot(plot_data_South, aes(order, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c","#5ec962")) +
  theme_classic()

ps
ggsave("Admixture_MiddleEastIncluded_South_K=4.pdf", plot = ps, width = 40, height = 10, unit = 'cm')



# MiddleEast
ME_Samples =Merged_origin[which(Merged_origin$OriginLessBroad=="Middle East"),1]
tbl_ME = tbl[which(tbl_ind[,1] %in% ME_Samples),]



plot_data_Me <- tbl_ME %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V4) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
                assingment_prob = max(prob)) %>% 
  dplyr::arrange(likely_assignment, desc(assingment_prob)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(id = forcats::fct_inorder(factor(id)))



plot_data_Me$order = rep("NA",nrow(plot_data_Me))
count=1
for (i in seq(1,length(levels(plot_data_Me$id)))){
  old_level = levels(plot_data_Me$id)[i]
  idx=which(plot_data_Me$id==old_level)
  plot_data_Me$order[idx]=count
  count=count+1
}
plot_data_Me$order=as.factor(as.integer(as.character(plot_data_Me$order)))


pme = ggplot(plot_data_Me, aes(order, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c","#5ec962")) +
  theme_classic()

pme
ggsave("Admixture_MiddleEastIncluded_MiddleEast_K=4.pdf", plot = pme, width = 40, height = 10, unit = 'cm')


# No location
NoLocation =Merged_origin[which(is.na(Merged_origin$OriginLessBroad=="South")),1]
tbl_NoLoc= tbl[which(tbl_ind[,1] %in% NoLocationIdx1),]

NoLocation_Samples = tbl[which(tbl_ind[,1] %in% NoLocation),]

plot_data_NoLoc <- NoLocation_Samples %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V4) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
                assingment_prob = max(prob)) %>% 
  dplyr::arrange(likely_assignment, desc(assingment_prob)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(id = forcats::fct_inorder(factor(id)))

plot_data_NoLoc$order = rep("NA",nrow(plot_data_NoLoc))
count=1
for (i in seq(1,length(levels(plot_data_NoLoc$id)))){
  old_level = levels(plot_data_NoLoc$id)[i]
  idx=which(plot_data_NoLoc$id==old_level)
  plot_data_NoLoc$order[idx]=count
  count=count+1
}
plot_data_NoLoc$order=as.factor(as.integer(as.character(plot_data_NoLoc$order)))


pnl = ggplot(plot_data_NoLoc, aes(order, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c","#5ec962")) +
  theme_classic()

pnl
ggsave("Admixture_MiddleEastIncluded_NoLocation_K=4.pdf", plot = pnl, width = 40, height = 10, unit = 'cm')





# bring together
Mixed_3 = rbind(plot_data_North,plot_data_South,plot_data_Asia,plot_data_Me,plot_data_NoLoc)
levels(Mixed_3$order)
Mixed_3$Location = rep("NA",nrow(Mixed_3))
Mixed_3$Location[1:560]="North"
Mixed_3$Location[561:1216]="South"
Mixed_3$Location[1217:1508]="Asia"
Mixed_3$Location[1509:1624]="Middle East"


Mixed_3$order = unfactor(Mixed_3$order)
Mixed_3$order[561:1216] = Mixed_3$order[561:1216] +1000 # larger level numbers for south
Mixed_3$order[1217:1508] = Mixed_3$order[1217:1508] +10000 # even larger level numbers for asia
Mixed_3$order[1509:1624] = Mixed_3$order[1509:1624] +100000 # even even larger level numbers for middle east
Mixed_3$order = as.factor(Mixed_3$order)

Mixed_3$Location = as.factor(Mixed_3$Location)
Mixed_3$Location <- factor(Mixed_3$Location, levels = c("North", "South", "Middle East","Asia","NA"))

K3AllGeo = ggplot(Mixed_3, aes(order, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c","#5ec962")) +
  theme_classic()  +facet_grid(. ~ Location,scales = "free", space = "free_x")  +  theme(axis.text.x = element_text(angle = 90))


ggsave("Admixture_ByGeographicLocation_K=4.pdf", plot = K3AllGeo, width = 100, height = 20, unit = 'cm')




####################################################
####################################################
######                K = 15                  ###### 
####################################################
####################################################

tbl=read.table("genotypesADMIXTURE_MAFfiltered.15.Q")


plot_data <- tbl %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V15) %>% 
  dplyr::group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
                assingment_prob = max(prob)) %>% 
  dplyr::arrange(likely_assignment, desc(assingment_prob)) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(id = forcats::fct_inorder(factor(id)))

p15 = ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col(width=1) + scale_fill_manual(values = c("#440154","#fde725","#21918c","#5ec962","orange","royalblue","lightslateblue","lightyellow","gray47","salmon1","olivedrab1","black","yellow3","skyblue1","midnightblue")) +
  theme_classic()

p15
ggsave("Admixture_Vir_K=15.pdf", plot = p15, width = 40, height = 10, unit = 'cm')


