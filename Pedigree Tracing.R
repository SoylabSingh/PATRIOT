#####For complex case (multiple pedigrees examined)#####
setwd('/Users/User/Desktop/Imputed')
markerfile <- read.csv('SoyNAMv2markers.csv')
Z <- read.csv('SoyNAMpedigrees.csv',header = T)

A <- markerfile[,1:3]
J <- as.data.frame(Z[,1])
for(x in 1:nrow(Z)){
  S <- as.data.frame(Z[x,])
  Self <- as.character(S[1,1])
  P1 <-as.character(S[1,2])
  P2 <- as.character(S[1,3])
  L <- markerfile[,which(colnames(markerfile) %in% c(Self,P1,P2))]
  H <- L[,c(Self,P1,P2)]
  rownames(H) <- A$ID
  for(i in 1:nrow(H)){
    Q = H[i,]
    b <- Q[,1] == Q[,2]
    c <- Q[,1] == Q[,3]
    d <- Q[,2] == Q[,3]
    P <- as.list(colnames(Q))
    E <- P[2]
    G <- P[3]
    if(b==TRUE & c==TRUE){H[i,4] = paste(E,G,sep=",",collapse="/")}
    else if(b==TRUE & c==FALSE){H[i,4] = E}
    else if(b==FALSE & c==TRUE){H[i,4] = G}
    else if(b==FALSE & c==FALSE){H[i,4] = 'Unknown'}
  }
  Y <- as.data.frame(H[,4])
  colnames(Y) <- Self
  A <- cbind(A,Y)
}
write.csv(A,"testbedtracedmarkers.csv")

#####For complex case (multiple pedigrees examined, NAM permutation)#####
setwd('/Users/User/Desktop/Imputed')
markerfile <- read.csv('SoyNAMv2markers.csv')
Z <- read.csv('SoyNAMpedigrees.csv',header = T)
A <- markerfile[,1:3]
J <- as.data.frame(Z[,1])
for(x in 1:nrow(Z)){
  S <- as.data.frame(Z[x,])
  Self <- as.character(S[1,1])
  P1 <-as.character(S[1,2])
  P2 <- as.character(S[1,3])
  L <- markerfile[,which(colnames(markerfile) %in% c(Self,P1,P2))]
  H <- L[,c(Self,P1,P2)]
  rownames(H) <- A$ID
  
  for(i in 1:nrow(H)){
    Q = H[i,]
    b <- Q[,1] == Q[,2]
    c <- Q[,1] == Q[,3]
    d <- Q[,2] == Q[,3]
    P <- as.list(colnames(Q))
    E <- P[2]
    G <- P[3]
    if(b==TRUE & c==TRUE){H[i,4] = paste(E,G,sep=",",collapse="/")}
    else if(b==TRUE & c==FALSE){H[i,4] = E}
    else if(b==FALSE & c==TRUE){H[i,4] = G}
    else if(b==FALSE & c==FALSE){H[i,4] = 'Unknown'}
  }
  Y <- as.data.frame(H[,4])
  colnames(Y) <- Self
  A <- cbind(A,Y)
}

write.csv(A,"pedigreestracedmarkersNAM09152020.csv")
###At this point, open the file and remove first column. I also renamed Chr and Pos columns.###
setwd('C:/Users/User/Desktop/Imputed/')
Markers_unimputed <- A
Markers_unimputed$CHROM <- as.factor(A$CHROM)
# identify number of levels
t_l=levels(Markers_unimputed$CHROM)
# convert all the variables in the column from factor to character
Markers_unimputed[] <- lapply(Markers_unimputed, as.character)

for(x in 1:length(t_l)){
  B <- t_l[x]
  C <- Markers_unimputed[which(Markers_unimputed$CHROM==B),]
  for(i in 4:ncol(C)){
    A <- unique(unlist(C[,i], use.names = FALSE))
    Line1 <- A[2]
    Line2 <- A[4]
    ParentList <- cbind.data.frame(Line1,Line2)
    BothPoss <- grep(',', A, value=TRUE)
    Unknown <- grep('Unknown', A, value=TRUE)
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
  write.csv(C,paste(B,".csv",sep=""))
}




#####For complex case (multiple pedigrees examined) SNP Correlation#####
setwd('/Users/jmshook/Desktop')
markerfile <- read.csv('PedMarkers.csv')
Z <- read.csv('PedTracetoDesktop.csv',header = T)
A <- markerfile[,1:3]
J <- as.data.frame(Z[,1])
for(x in 80:nrow(Z)){
  S <- as.data.frame(Z[x,])
  Self <- as.character(S[1,1])
  P1 <-as.character(S[1,2])
  P2 <- as.character(S[1,3])
  L <- markerfile[,which(colnames(markerfile) %in% c(Self,P1,P2))]
  H <- L[,c(Self,P1,P2)]
  rownames(H) <- A$ID
  
  for(i in 1:nrow(H)){
    Q = H[i,]
    b <- Q[,1] == Q[,2]
    c <- Q[,1] == Q[,3]
    d <- Q[,2] == Q[,3]
    P <- as.list(colnames(Q))
    E <- P[2]
    G <- P[3]
    if(b==TRUE & c==TRUE){H[i,4] = 'NA'}
    else if(b==TRUE & c==FALSE){H[i,4] = -1}
    else if(b==FALSE & c==TRUE){H[i,4] = 1}
    else if(b==FALSE & c==FALSE){H[i,4] = 'NA'}
  }
  Y <- as.data.frame(H[,4])
  colnames(Y) <- Self
  A <- cbind(A,Y)
}

write.csv(A,"pedigreestracedmarkers.csv")
#####SNP correlations peds
rownames(A) <- A$ID
install.packages('Hmisc')
library(Hmisc)
AB <- A[A=="NA"] <- NA
SNPCorr <- rcorr(as.matrix(t(A[,4:ncol(A)])),type="pearson")
gq <- as.matrix(SNPCorr[["r"]])
gq[1:5,1:5]
ZZTop <- data.frame(row=rownames(gq)[row(gq)[upper.tri(gq)]], 
           col=colnames(gq)[col(gq)[upper.tri(gq)]], 
           corr=gq[upper.tri(gq)])
write.csv(ZZTop,'SNPcorrspedtracing.csv')
gq2 <- as.data.frame(gq)
gq3 <- data.frame(lapply(gq2,function(x){gsub("NA","",x)}))
gq3 <- gq2[,which(colnames(gq2) %in% Boop$ID)]
gq3[-0.8 < gq3] <- NA
gq6 <- gq5[gq5 < 0.8] <- NA
