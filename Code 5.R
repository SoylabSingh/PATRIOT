#####Start new code#####
setwd('/Users/User/Desktop/PATRIOT')
#Read in imputed file
Finished4 <- read.csv('imputeexample.csv')
#Create empty df to handle everything after the loop finishes.
Finished7 <- data.frame()
#Create empty df to handle everything inside the loop
Finished6 <- data.frame()

#Skip first three columns (Chr, Pos, MarkerName) and start with first imputed progeny.
for(x in 4:ncol(Finished4)){
#Create object with column names as values
  col <- colnames(Finished4)
#Reset index to 1
  start1=1
#Loop through chromosomes
  for(y in 1:20){
#Subset imputed df to the chromosome under consideration
    Finished5 <- Finished4[which(Finished4$CHROM == y),]
#Reset index to 1 to begin new loop
    start1=1
#Loop through subsequent pairs of markers
    for(i in 2:nrow(Finished5)){
#If subsequent pairs of markers do not match (not inherited from same source)
      if(Finished5[i,x] != Finished5[i-1,x]){
#First column gets the name of the progeny being analyzed
        Finished6[1,1] <- col[x]
#Second column gets the current chromosome
        Finished6[1,2] <- Finished5$CHROM[start1]
#Third column gets the position of the first marker in the string of markers inherited from the same source
        Finished6[1,3] <- Finished5$POS[start1]
#Fourth column gets the position of the last marker to retain the initial ancestral source
        Finished6[1,4] <- Finished5$POS[i-1]
#Fifth column gets the name of the ancestral source
        Finished6[1,5] <- Finished5[i-1,x]
#Set the index to the first marker which was not inherited from the previous ancestral source
        start1 = i
#Add the row to the permanently saved df for writing later.         
        Finished7 <- rbind(Finished7,Finished6,stringsAsFactors=FALSE)
#Regenerate a clean df to work in
        Finished6 <- data.frame(stringsAsFactors = FALSE)}
#If final marker is reached
      if(i == nrow(Finished5)){
#First column gets the name of the progeny being analyzed
        Finished6[1,1] <- col[x]
#Second column gets the current chromosome
        Finished6[1,2] <- Finished5$CHROM[start1]
#Third column gets the position of the first marker in the string of markers inherited from the same source
        Finished6[1,3] <- Finished5$POS[start1]
#Fourth column gets the position of the last marker. This loop differs from the previous in order to capture the correct marker position for the last marker on the chromosome.
        Finished6[1,4] <- Finished5$POS[i]
#Fifth column gets the name of the ancestral source.        
        Finished6[1,5] <- Finished5[i,x]
#Add the row to the permanently saved df for writing later.         
        Finished7 <- rbind(Finished7,Finished6,stringsAsFactors=FALSE)
#Regenerate a clean df to work in        
        Finished6 <- data.frame(stringsAsFactors = FALSE)}}}}
#Write the finished df to file
write.csv(Finished7,"recombinationexample.csv")
