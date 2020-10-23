setwd('/Users/User/Desktop/Imputed')
markerfile <- read.csv('SoyNAMv2markers.csv')
#####LENGTH 42.449 BASED ON V1. V2 WOULD BE 42.080####
Markers_unimputed <- markerfile[1:4289,]
# identify number of levels
t_l=levels(Markers_unimputed$CHROM)
# convert all the variables in the column from factor to character
Markers_unimputed[] <- lapply(Markers_unimputed, as.character)
Finished <- data.frame()
######start forward loop#####
for(x in 1:20){
  C <- Markers_unimputed[which(Markers_unimputed$CHROM==x),]
  for(i in 4:ncol(C)){
    D <- unique(unlist(C[,i], use.names = FALSE))
    BothPoss <- grep(',', D, value=TRUE)
    Unknown <- grep('Unknown', D, value=TRUE)
    Impute <- list(BothPoss,Unknown)
    ParentList=unlist(strsplit(BothPoss,','))
    for(n in 1:(nrow(C))){
      if (C[n,i] %in% Impute){
        for (d in n:nrow(C))  {
          if (C[d,i] %in% ParentList)  {
            bottomparent <- C[d,i]
            break}
          bottomparent <- C[d,i]}
        for (e in n:1)  {
          if (C[e,i] %in% ParentList)  {
            topparent <- C[e,i]
            break}
          topparent <- C[e,i]}
        if (bottomparent == topparent){
          C[n,i] = topparent}
        else if(bottomparent != topparent){
          C[n,i]=C[n,i]}
        else if(C[n,i] %in% ParentList){ C[n,i] <- C[n,i]}}
    }
  }
  j <- C
  Finished <- rbind(Finished,j)
  rm(bottomparent,bottomparent2,topparent,ParentList,Impute,BothPoss,D)
}
write.csv(Finished,'NAMmiddleimpute09152020.csv')
View(Finished)
#####End middle loop#####
#####Start forward loop#####
setwd('C:/Users/User/Desktop/Imputed')
Finished <- read.csv('NAMmiddleimpute09152020.csv')
Finished2 <- data.frame()
for (x in 1:20) {
  C <- Finished[which(Finished$CHROM==x),]
  for (i in 4:ncol(C)) {
    D <- levels(C[,i])
    BothPoss <- grep(',',D,value = TRUE)
    Unknown <- grep('Unknown',D,value = TRUE)
    Impute <- list(BothPoss,Unknown)
    ParentList = unlist(strsplit(BothPoss,','))
    if (C[1,i] %in% Impute) {
      for (q in 1:nrow(C)) {
        if (C[q,i]%in% ParentList) {
          bottomparent2 <- C[q,i]
        break}
        bottomparent2 <- C[q,i]
      }    
      for (g in q:1) {
        C[g,i] <- bottomparent2}
    }
    else if(C[1,i] %in% ParentList){
      bottomparent2 <- C[1,i]        
    }
    }
  j <- C
  Finished2 <- rbind(Finished2,j)
}
write.csv(Finished2,"NAMtopimpute09152020.csv")#no errors so far
#####Start backward loop#####
beta <- as.data.frame(Finished2)
beta$CHROM <- as.numeric(beta$CHROM)
beta$POS <- as.numeric(beta$POS)
Finished3 <- beta[order(beta$CHROM,-beta$POS),]
write.csv(Finished3,"ggg09152020.csv")
Finished3 <- read.csv("ggg09152020.csv")
Finished3$X <- NULL
Finished4 <- data.frame()
for (x in 1:20) {
  C <- Finished3[which(Finished3$CHROM==x),]
  for (i in 4:ncol(C)) {
    D <- levels(C[,i])
    BothPoss <- grep(',',D,value = TRUE)
    Unknown <- grep('Unknown',D,value = TRUE)
    Impute <- list(BothPoss,Unknown)
    ParentList = unlist(strsplit(BothPoss,','))
    if (C[1,i] %in% Impute) {
      for (q in 1:nrow(C)) {
        if (C[q,i]%in% ParentList) {
          bottomparent2 <- C[q,i]
          break}
        bottomparent2 <- C[q,i]
      }    
      for (g in q:1) {
        C[g,i] <- bottomparent2}
    }
    else if(C[1,i] %in% ParentList){
      bottomparent2 <- C[1,i]        
    }
  }
  j <- C
  Finished4 <- rbind(Finished4,j)
}
#####End backward loop#####
write.csv(Finished4,"tracedNAMnew09152020.csv")