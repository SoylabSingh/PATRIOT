#####For complex case (multiple pedigrees examined)#####
setwd('/Users/User/Desktop/PATRIOT')
#Read marker file. May be stored as '0,1,2', 'A,B,H', or 'A,T,C,G,H'#
markerfile <- read.csv('SoyNAMv2markers.csv')
#Read in Pedigree file. Progeny name in first column, female parent in second, male parent in third.#
Z <- read.csv('SoyNAMpedigrees.csv',header = T)
#Subset the Chromosome, Position, and Marker Name columns#
A <- markerfile[,1:3]
#Create list of progeny to run the pipeline on.#
J <- as.data.frame(Z[,1])
#Loop through each progeny.#
for(x in 1:nrow(Z)){
#Comment the next line to run all, or comment the preceding line to run all progeny.#
#for(x in 1:5){
#Get pedigree for xth progeny.#
  S <- as.data.frame(Z[x,])
#Extract progeny name.# 
  Self <- as.character(S[1,1])
#Extract female parent name.#  
  P1 <-as.character(S[1,2])
#Extract male parent name.#
  P2 <- as.character(S[1,3])
#Subset marker files to those in the immediate pedigree.#
  L <- markerfile[,which(colnames(markerfile) %in% c(Self,P1,P2))]
#Order columns.#
  H <- L[,c(Self,P1,P2)]
#Set marker name as rowname.#
  rownames(H) <- A$ID
#Loop through markers.#
  for(i in 1:nrow(H)){
#Create subset of only the ith position within the immediate pedigree.#
    Q = H[i,]
#Does the progeny match the female parent at this position?#
    b <- Q[,1] == Q[,2]
#Does the progeny match the male parent at this position?#    
    c <- Q[,1] == Q[,3]
#Are the parents fixed at this position?#
    d <- Q[,2] == Q[,3]
#List progeny and parents#
    P <- as.list(colnames(Q))
#Female parent name.#
    E <- P[2]
#Male parent name.#
    G <- P[3]
#If marker is fixed for both parents and the progeny, replace the marker with 'Parent1/Parent2'.#
    if(b==TRUE & c==TRUE){H[i,4] = paste(E,G,sep=",",collapse="/")}
#If marker matches only female parent, replace the marker with 'Parent1'.#
    else if(b==TRUE & c==FALSE){H[i,4] = E}
#If marker matches only male parent, replace the marker with 'Parent2'.#
    else if(b==FALSE & c==TRUE){H[i,4] = G}
#Otherwise, write that the marker source is 'Parent1/Parent2'.#    
    else if(b==FALSE & c==FALSE){H[i,4] = paste(E,G,sep=",",collapse="/")}
  }
#Write the new traced marker names to a dataframe.#
  Y <- as.data.frame(H[,4])
#Set column names to progeny name.#
  colnames(Y) <- Self
#Column bind to 1) column, position, marker name df and 2) previous traced markers.#
  A <- cbind(A,Y)
}
#Write the final file for further use.#
write.csv(A,"testbedtracedmarkers.csv")
#Free up space for further analyses.#
rm(b,c,d,E,G,H,i,J,L,markerfile,P,P1,P2,Q,S,Self,x,Y,Z)
#####Author's note#####
#Working in an inbred crop, heterozygous markers are not rigorously handled in this script. 
#Modifications will be needed for adoption with heterozygous species, or to handle phased chromosomes.#