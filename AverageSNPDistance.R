
# A script to calculate the distances between two adjacent SNPs


vcf_file = read.table("Mapped_HighQualSNPs2022-01-06_imputed.vcf",stringsAsFactors=F)

positions = vcf_file[,1:2]

n = nrow(positions)
cat(paste("Found",n,"SNPs",sep=" "))


# chromosome-wise dataframes
Chr1S = positions[which(positions[,1]=="chr1S"),]
Chr1L = positions[which(positions[,1]=="chr1L"),]
Chr2 = positions[which(positions[,1]=="chr2"),]
Chr3 = positions[which(positions[,1]=="chr3"),]
Chr4 = positions[which(positions[,1]=="chr4"),]
Chr5 = positions[which(positions[,1]=="chr5"),]
Chr6 = positions[which(positions[,1]=="chr6"),]



Distance_AdSNPs <- function(df){
  #function that calculates distance between all adjacent SNPs and reports 
  #maximum, minimum and average as well
  
  # make df with n rows
  m <- matrix(0, ncol = 1, nrow = (nrow(df)-1))
  m <- data.frame(m)
  colnames(m) = "distance_between_two_adjacent_SNPs"
  
  df[,2] = as.numeric(as.character(df[,2]))
  df_ordered = df[order(df[,2]),]
  
  for (i in seq(1:((nrow(df)-1)))){
    
    cat(paste("At i =",i,"out of", n, "SNPs", sep=" "))
    distance <- df_ordered[i+1,2]-df_ordered[i,2]
    m$distance_between_two_adjacent_SNPs[i] <- distance
  }
  
  average_dist = mean(m$distance_between_two_adjacent_SNPs)
  max_dist = max(m$distance_between_two_adjacent_SNPs)
  min_dist = min(m$distance_between_two_adjacent_SNPs)
  
  new_m = rbind(m,average_dist,max_dist,min_dist)  
  return(new_m)
}


Matrix = Distance_AdSNPs(Chr1S)
write.table(Matrix,"Chr1S_distances.txt",sep="\t",row.names=F,col.names=F,quote=F)

Matrix = Distance_AdSNPs(Chr1L)
write.table(Matrix,"Chr1L_distances.txt",sep="\t",row.names=F,col.names=F,quote=F)

Matrix = Distance_AdSNPs(Chr2)
write.table(Matrix,"Chr2_distances.txt",sep="\t",row.names=F,col.names=F,quote=F)

Matrix = Distance_AdSNPs(Chr3)
write.table(Matrix,"Chr3_distances.txt",sep="\t",row.names=F,col.names=F,quote=F)

Matrix = Distance_AdSNPs(Chr4)
write.table(Matrix,"Chr4_distances.txt",sep="\t",row.names=F,col.names=F,quote=F)

Matrix = Distance_AdSNPs(Chr5)
write.table(Matrix,"Chr5_distances.txt",sep="\t",row.names=F,col.names=F,quote=F)

Matrix = Distance_AdSNPs(Chr6)
write.table(Matrix,"Chr6_distances.txt",sep="\t",row.names=F,col.names=F,quote=F)



