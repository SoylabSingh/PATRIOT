#Bring in 'A' dataframe from Pedigree Tracing script.#
#If you need to, read it it from file.#
#setwd('/Users/User/Desktop/PATRIOT')
#A <- read.csv('testbedtracedmarkers.csv')
#Remove the first column (numeric rownames from saving).#
#A$X <- NULL

#Identify number of levels (chromosomes).#
t_l=levels(A$CHROM)
#Convert all the variables in the column from factor to character.#
A[] <- lapply(A, as.character)
#Create new dataframe to fill in.#
Finished <- data.frame()
#Pick one chromosome at a time.#
#With work, this could be parallel-ized.#
for(x in 1:20){
  C <- A[which(A$CHROM==x),]
#Loop through progeny.#  
  for(i in 4:ncol(C)){
#Get list of all factor classes.#    
    D <- unique(unlist(C[,i], use.names = FALSE))
#Identify 'Parent1/Parent2' type unknowns for imputing.#
    BothPoss <- grep(',', D, value=TRUE)
#Set a class for 'Unknown's.#
    Unknown <- grep('Unknown', D, value=TRUE)
#Create list of factors to impute (Unknown and P1/P2).#
    Impute <- list(BothPoss,Unknown)
#Parents can be from either of the halves of Bothposs. Split them at the comma.#
    ParentList=unlist(strsplit(BothPoss,','))
#Record class of first (Red) and last (Blue) marker on current chromosome for ith line.
    Red <- C[1,i]
    Blue <- C[nrow(C),i]
#If first marker in Impute
    if(Red %in% Impute){
#Loop through following markers      
      for(p in 2:nrow(C)){
#Until you reach a marker that is known to come from one of the parents
        if(C[p,i] %in% ParentList){
#Save the name of that parent.#
          bottomparent2 <- C[p,i]
#Loop through the first p markers          
          for(o in 1:p){
#And replace with that parent's name.#            
            C[o,i] <- bottomparent2}
#Make sure not to loop through additional markers after finding the first parental occurrence.#          
         break}}}
#If last marker in Impute
    if(Blue %in% Impute){
#Loop through preceding markers      
      for(w in nrow(C):1){
#Until you reach a marker that is known to come from one of the parents
        if(C[w,i] %in% ParentList){
#Save the name of that parent.#
          topparent2 <- C[w,i]
#Loop through the markers following marker w  
          for(e in w:nrow(C)){
#And replace with that parent's name.#  
            C[e,i] <- topparent2}
#Make sure not to loop through additional markers after finding the first parental occurrence.#  
        break}}}
    
#Loop through markers.#
    for(n in 1:(nrow(C))){
#If the nth marker for progeny i is either Unknown (doesn't match either parent) or fixed (P1/P2 style format, caused by both parents as well as the progeny fixed for the same allele)... 
      if (C[n,i] %in% Impute){
#Start iterating downward...
        for (d in n:nrow(C)){
#until you reach the next marker with known parent.#
          if (C[d,i] %in% ParentList)  {
#Set that marker's class as the "bottom parent".
#It's not a great name, but hey, it's closer to the bottom of the dataframe.#
            bottomparent <- C[d,i]
            break}
#When testing, this next line fixed issues I was having, but I do NOT understand why.#
          bottomparent <- C[d,i]}
#Start iterating upwards from the marker to be imputed...
        for (e in n:1){
#until you reach the next marker with known parent.
          if (C[e,i] %in% ParentList){
#Set that marker's class as the "top parent".#
#Suggestions for better naming for top/bottom parent, as well as any other variables is welcome.#
            topparent <- C[e,i]
            break}
#When testing, this next line fixed issues I was having, but I do NOT understand why.#
          topparent <- C[e,i]}
#If the top and bottom parents are the same...
        if (bottomparent == topparent){
#Replace the cell to be imputed with their name. I chose topparent, but bottomparent here would have the same effect.#
          C[n,i] = topparent}
#If the top and bottom parents are different...
        else if(bottomparent != topparent){
#Don't do anything.#
          C[n,i]=C[n,i]}
#If the nth marker for the ith progeny was already known, keep it.#
        else if(C[n,i] %in% ParentList){ C[n,i] <- C[n,i]}}}}
#After a chromosome is finished, rename and append it to the Finished dataframe.
  j <- C
#Append it to said dataframe.#
  Finished <- rbind(Finished,j)
#Clear these out to prevent them from accidentally carrying over to next progeny.#
  rm(bottomparent,topparent,ParentList,Impute,BothPoss,D,topparent2,bottomparent2)}
#Write to file if desired. This can serve as a sanity check when coupled with saving the initial (Tracing) file.#
write.csv(Finished,'imputeexample.csv')
#####Author's notes#####
#The first 10 markers of Chromosome 1 in the example demonstrate the imputation rules for 
#1) multiple consecutive markers (2-5) flanked by the same parent, 
#2) a single marker flanked by the same parent (8), and 
#3) a chunk of markers (or single marker) flanked by different parents.#
#Chromosome 3 of the third progeny (DS11.02006) demonstrates the rules for end-of-chromosome imputation when the first or last marker on the chromosome are not mapped to exactly one parent.
#The first (last) marker on the chromosome which is mapped to exactly one parent is used to impute the preceding (following) markers.#