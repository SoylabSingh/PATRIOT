#Install packages needed
install.packages("tidyr")
install.packages("stringi")
install.packages("stringr")
install.packages("plyr")
install.packages("dplyr")
install.packages('rrBLUP')
install.packages('data.table')
library(tidyr)
library(stringi)
library(data.table)
library(rrBLUP)
library(plyr)
library(dplyr)
library(stringr)

#####Loop through GP#####
#####Initiate the environment with marker files#####
setwd('/Users/User/Desktop/PATRIOT')
A <- read.csv('imputeexample.csv')
A2 <- read.csv('SoyNAMv2markersnumeric.csv')
#Read in Pedigree file. Progeny name in first column, female parent in second, male parent in third.#
Z <- read.csv('SoyNAMpedigrees.csv',header = T)
#Create list of all possible parents from pedigree#
Unique <- unique(c(Z[,2],Z[,3]))
#Create list of possible parent combinations from pedigree#
Predup <- unique(Z[,2:3])
#Coerce parental combinations to same format as Pedigree Tracing.R#
Dupe <- vector(mode="list",length=nrow(Predup))
for(j in 1:nrow(Predup)){
  Dupe[j] <- paste(Predup[j,1],",",Predup[j,2])}
Dupe <- gsub(" ","",Dupe)
#Set rownames
rownames(A) <- A$SNP
#Remove Chr, Pos, MarkerName columns from 'traced' file
A$SNP <- NULL
A$CHROM <- NULL
A$POS <- NULL
#Transpose for use in calculating allele effect estimates 
A1 <- data.frame(t(A))
#Set rownames
rownames(A2) <- A2$SNP
#Remove Chr, Pos, MarkerName columns from raw marker file (0,1,2 format)
A2$SNP <- NULL
A2$CHROM <- NULL
A2$POS <- NULL
#Transpose
A3 <- data.frame(t(A2))
#Read in phenotypic records- in this case, multi-trait, multi-environment.
df <- read.csv('NAMFullPheno.csv')
#Generate list of environments
enviros=unique(df$Env)
#Clear space
rm(A2,A)
#Remove NAM parents from consideration in GS model
A3 <- A3[41:5189,] 
#Initiate df to store correlations in
Recordcorrelations <- data.frame()
#Garbage collection
gc()
gc()
gc()

#We care about A1,A3,df,enviros. Others are superfluous.
#Loop through environments if multi-environment data.
for(q in 1:length(enviros)){
#Get environment name
  p <- enviros[q]
#Subset phenotypic records to those from the qth environment
  s <- df[which(df$Env==p),]
#Set rownames to line names
  rownames(s) <- s$CorrectedStrain
#Remove line name column
  s$CorrectedStrain <- NULL
######Change based on trait#####  
  s <- s %>% drop_na(Yldkgha) 
#Subset lines with phenotype records to those in traced marker file  
  df1 <- s[which(rownames(s) %in% rownames(A1)),]
#Subset traced marker file to those in phenotyic records file
  A2 <- A1[which(rownames(A1) %in% rownames(df1)),]
#Subset lines with phenotype records to those in numeric marker file
  df2 <- s[which(rownames(s) %in% rownames(A3)),]
#Subset numeric marker file to those with phenotypic records
  A4 <- A3[which(rownames(A3) %in% rownames(df2)),]
#Determine number of lines with both phenotypic and genotypic records
  u <- nrow(df1)
#Subset out the trait of interest (only necessary for multi-trait records)
  Pheno <- df1 %>% select(Yldkgha) #change based on trait
#For evaluation only, we subset 80% of lines to use as training
  training_entries <- as.matrix(sample(1:u,floor(u*0.8)))
#The lines not used for training are used as testing (unseen). 
#This step can be modified for production use to take those with phenotypic records as training and those with only genotypic records as testing.
  testing_entries <- setdiff(1:u,training_entries)
#Subset phenotypic records for training set  
  Pheno_training_data = as.matrix(Pheno[training_entries,])
#Set the column name for the phenotypic records so that it will be consistent for downstream steps.
  colnames(Pheno_training_data) <- "aeematrix"
#Subset phenotypic records for testing set
  Pheno_testing_data = as.matrix(Pheno[testing_entries,])
#Subset marker data for traced training set
  SNP_training_data = as.matrix(A2[training_entries,])
#Subset marker data for traced testing set
  SNP_testing_data = as.matrix(A2[testing_entries,])  
#Subset marker data for numeric training set
  SNP_training_data2 = as.matrix(A4[training_entries,])
#Subset marker data for numeric testing set
  SNP_testing_data2 = as.matrix(A4[testing_entries,])
#Attach phenotypic data to the end of the traced training marker set
  B <- cbind(SNP_training_data,as.data.frame(Pheno_training_data))
#Create df to store allele effect estimates
  aeematrix <- data.frame()
#Loop through each marker
  for(a in 1:(ncol(B)-1)) {
#Loop through known parents
    for (h in 1:length(Unique)) {
#Calculate the mean phenotype of lines with current allele coming from specific parent compared to entire population    
  aeematrix[h,a] <- mean(B[,ncol(B)][which(B[,a]==Unique[h])])-mean(B[,ncol(B)])}}
  rownames(aeematrix) <- Unique
  colnames(aeematrix) <- colnames(B[1:(ncol(B)-1)])
#NA handling
  aeematrix[aeematrix=="NaN"] <- 0
  PostParents <- nrow(aeematrix)
#Loop through markers#
  for(b in 1:(ncol(B)-1)) {
#Loop through parental combinations#
    for (c in 1:length(Dupe)) {
#Get the first potential list of parent pairs#
      TwoParents <- unlist(strsplit(Dupe[[c]],','))
#Extract the first parent name#
      FirstPar <- TwoParents[[1]]
#Extract the second parent name#
      SecondPar <- TwoParents[[2]]
#Record the AEE of the first parent
      rose <- aeematrix[FirstPar,b]
#Record the AEE of the second parent      
      rose2 <- aeematrix[SecondPar,b]
#Average the two and save as new column in aeematrix
      aeematrix[c+PostParents,b] <- 0.5*(rose+rose2)
}}
FullRowNameList <- c(Unique,Dupe)  
rownames(aeematrix) <- FullRowNameList
for(Zed in 1:ncol(aeematrix)){
#Create df to look up AEEs and replace traced markers with
    replacewith <- as.data.frame(aeematrix[,a])
#Add the name associated with each class. It is critical to have the names be in the proper order.
    replacewith$Lookup <- rownames(aeematrix)
#Set column name
    colnames(replacewith) <- c("AEE","Lookup")   
    replacement <- replacewith[seq(dim(replacewith)[1],1),]
#Replace traced marker classes with AEEs in training dataframe
    SNP_training_data[,Zed] <- stri_replace_all_fixed(SNP_training_data[,Zed], pattern = replacement$Lookup, replacement = replacement$AEE, vectorize_all = FALSE)
#Replace traced marker classes with AEEs in testing dataframe
    SNP_testing_data[,Zed] <- stri_replace_all_fixed(SNP_testing_data[,Zed], pattern = replacement$Lookup, replacement = replacement$AEE, vectorize_all = FALSE)
}
#Convert to data.frame  
  SNP_testing_data <- data.frame(SNP_testing_data)
#Convert to data.frame
  SNP_training_data <- data.frame(SNP_training_data)
#Next four lines- write and read to get AEE data in a format suitable for rrBLUP  
  fwrite(SNP_training_data,"train.temp")
  b1 <- fread("train.temp",colClasses = "numeric")
  fwrite(SNP_testing_data,"test.temp")
  b2 <- fread("test.temp",colClasses = "numeric")
#Train a model with the AEEs as a replacement for raw marker data
  trained_model <- mixed.solve(y=Pheno_training_data,Z=b1)
#Record marker effects
  marker_effects <- as.matrix(trained_model$u)
#Generate BLUEs from the trained model
  BLUE <- as.vector(trained_model$beta)
#Calculate sum of marker effects for training set
  predicted_train = as.matrix(b1) %*% marker_effects
#Calculate sum of marker effects for testing set  
  predicted_test = as.matrix(b2) %*% marker_effects
#Add in the BLUEs to generate predictions for the training set  
  predicted_train_result <- as.vector((predicted_train[,1])+BLUE)
#Add in the BLUEs to generate predictions for the testing set
  predicted_test_result <- as.vector((predicted_test[,1])+BLUE)
#Do the same as above with the raw marker data (0,1,2 format)  
  trained_model2 <- mixed.solve(y=Pheno_training_data,Z=SNP_training_data2)
  marker_effects2 <- as.matrix(trained_model2$u)
  BLUE2 <- as.vector(trained_model2$beta)
  predicted_train2 = as.matrix(SNP_training_data2) %*% marker_effects2
  predicted_test2 = as.matrix(SNP_testing_data2) %*% marker_effects2
  predicted_train_result2 <- as.vector((predicted_train2[,1])+BLUE2)
  predicted_test_result2 <- as.vector((predicted_test2[,1])+BLUE2)
#First column gets the environment name  
  Recordcorrelations[q,1] <- p
#Second column gets the correlation between predicted and observed values for PATRIOT-rrBLUP testing set
  Recordcorrelations[q,2] <- cor(as.vector(Pheno_testing_data),predicted_test_result,use="complete")
#Third column gets the correlation between predicted and observed values for PATRIOT-rrBLUP training set
  Recordcorrelations[q,3] <- cor(as.vector(Pheno_training_data),predicted_train_result,use="complete")
#Fourth column gets the correlation between predicted and observed values for standard rrBLUP testing set 
  Recordcorrelations[q,4] <- cor(as.vector(Pheno_testing_data),predicted_test_result2,use="complete")
#Fifth column gets the correlation between predicted and observed values for standard rrBLUP training set
  Recordcorrelations[q,5] <- cor(as.vector(Pheno_training_data),predicted_train_result2,use="complete")
#Optionally, save AEE matrix.
  write.csv(aeematrix,"aeematrix1.csv")
  }
#Label columns to keep track of which column refers to which correlation
colnames(Recordcorrelations) <- c("Environment","TestingAEE","TrainingAEE","TestingStandard","TrainingStandard")
#Write to file
write.csv(Recordcorrelations,'Standard vs AEEYld21.csv')