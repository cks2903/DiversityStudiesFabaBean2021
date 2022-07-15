library(UpSetR)
library(ggplot2)


input <- c(
  Ohana = 1562,
  pcadapt = 304,
  bayescan = 12,
  "Ohana&pcadapt" = 29,
  "Ohana&bayescan" = 0,
  "pcadapt&bayescan" = 1,
  "Ohana&pcadapt&bayescan" = 5
)
pdf(file="/Volumes/NAT_MBG-PMg/Cathrine/Profaba/Meta_analysis_DiversityPanel_Jan2022/The_Diversity_Panel/20220314_SegregatingMarkersUnderSelection/OutlierDetectionUpsetplot.pdf",width = 20,height=6.5) # or other device
p <- upset(fromExpression(input),
     nintersects = 50, 
     nsets = 10, 
     order.by = "freq", 
     decreasing = T, 
     mb.ratio = c(0.6, 0.4),
     number.angles = 0, 
     text.scale = 2.8, 
     point.size = 2.8, 
     line.size = 1
)
p
dev.off()


