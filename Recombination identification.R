#####Start new code#####
setwd('/Users/User/Desktop/Imputed')
Finished4 <- read.csv('allimputed2.csv')
Finished7 <- data.frame()
Finished6 <- data.frame()
for(x in 4:ncol(Finished4)){
  col <- colnames(Finished4)
  start1=1
  for(y in 1:20){
    Finished5 <- Finished4[which(Finished4$CHROM == y),]
    start1=1
    for( i in 2: nrow(Finished5)){
      if(Finished5[i,x] != Finished5[i-1,x]){
        Finished6[1,1] <- col[x]
        Finished6[1,2] <- Finished5$CHROM[start1]
        Finished6[1,3] <- Finished5$POS[start1]
        Finished6[1,4] <- Finished5$POS[i-1]
        Finished6[1,5] <- Finished5[i-1,x]
        start1 = i
        Finished7 <- rbind(Finished7,Finished6,stringsAsFactors=FALSE)
        Finished6 <- data.frame(stringsAsFactors = FALSE)}
      if(i == nrow(Finished5)){
        Finished6[1,1] <- col[x]
        Finished6[1,2] <- Finished5$CHROM[start1]
        Finished6[1,3] <- Finished5$POS[start1]
        Finished6[1,4] <- Finished5$POS[i]
        Finished6[1,5] <- Finished5[i,x]
        Finished7 <- rbind(Finished7,Finished6,stringsAsFactors=FALSE)
        Finished6 <- data.frame(stringsAsFactors = FALSE)}
    }
  }
}
write.csv(Finished7,"recombinationfullpanel.csv")
