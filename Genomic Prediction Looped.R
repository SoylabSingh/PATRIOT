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
#####Initiate the environment with needed files#####
setwd('/Users/User/Desktop/Imputed')
A <- read.csv('tracedNAMnew02162020.csv')
A2 <- read.csv('SoyNAMv2markersnumeric.csv')
rownames(A) <- A$SNP
A$SNP <- NULL
A$CHROM <- NULL
A$POS <- NULL
A1 <- data.frame(t(A))
rownames(A2) <- A2$SNP
A2$SNP <- NULL
A2$CHROM <- NULL
A2$POS <- NULL
A3 <- data.frame(t(A2))
df <- read.csv('NAMFullPheno.csv')
enviros=levels(df$Env)
rm(A2,A)
A3 <- A3[41:5189,] 
Recordcorrelations <- data.frame()
gc()
gc()
gc()

#We care about A1,A3,df,enviros. Others are superfluous.
for(q in 1:length(enviros)){
  p <- enviros[q]  
  s <- df[which(df$Env==p),]
  rownames(s) <- s$CorrectedStrain
  s$CorrectedStrain <- NULL
  s <- s %>% drop_na(DaystoMat) #change based on trait
  df1 <- s[which(rownames(s) %in% rownames(A1)),]
  A2 <- A1[which(rownames(A1) %in% rownames(df1)),]
  df2 <- s[which(rownames(s) %in% rownames(A3)),]
  A4 <- A3[which(rownames(A3) %in% rownames(df2)),]
  u <- nrow(df1)
  Pheno <- df1 %>% select(DaystoMat) #change based on trait
  training_entries <- as.matrix(sample(1:u,floor(u*0.8)))
  testing_entries <- setdiff(1:u,training_entries)
  Pheno_training_data = as.matrix(Pheno[training_entries,])
  colnames(Pheno_training_data) <- "aeematrix"
  Pheno_testing_data = as.matrix(Pheno[testing_entries,])
  SNP_training_data = as.matrix(A2[training_entries,])
  SNP_testing_data = as.matrix(A2[testing_entries,])  
  SNP_training_data2 = as.matrix(A4[training_entries,])
  SNP_testing_data2 = as.matrix(A4[testing_entries,])
  B <- cbind(SNP_training_data,as.data.frame(Pheno_training_data))
  aeematrix <- data.frame()
  for (a in 1:(ncol(B)-1)) {
    dfIA3023 <- B[which(B[,a]=="IA3023"),]
    dfCL0J09546 <- B[which(B[,a]=="CL0J095.4.6"),]
    dfHS63976 <- B[which(B[,a]=="HS6.3976"),]
    dfLD015907 <- B[which(B[,a]=="LD01.5907"),]
    dfLD024485 <- B[which(B[,a]=="LD02.4485"),]
    dfLD029050 <- B[which(B[,a]=="LD02.9050"),]
    dfLG032979 <- B[which(B[,a]=="LG03.2979"),]
    dfLG033191 <- B[which(B[,a]=="LG03.3191"),]
    dfLG044717 <- B[which(B[,a]=="LG04.4717"),]
    dfLG054292 <- B[which(B[,a]=="LG05.4292"),]
    dfLG054317 <- B[which(B[,a]=="LG05.4317"),]
    dfLG054464 <- B[which(B[,a]=="LG05.4464"),]
    dfLG054832 <- B[which(B[,a]=="LG05.4832"),]
    dfNE3001 <- B[which(B[,a]=="NE3001"),]
    dfPI398881 <- B[which(B[,a]=="PI398881"),]
    dfPI404188A <- B[which(B[,a]=="PI404188A"),]
    dfPI427136 <- B[which(B[,a]=="PI427136"),]
    dfPI437169B <- B[which(B[,a]=="PI437169B"),]
    dfPI518751 <- B[which(B[,a]=="PI518751"),]
    dfPI561370 <- B[which(B[,a]=="PI561370"),]
    dfPI574486 <- B[which(B[,a]=="PI574486"),]
    dfPI595362 <- B[which(B[,a]=="PI595362"),]
    dfPI598124 <- B[which(B[,a]=="PI598124"),]
    dfPI602995 <- B[which(B[,a]=="PI602995"),]
    dfPI615553 <- B[which(B[,a]=="PI615553"),]
    dfPI633730 <- B[which(B[,a]=="PI633730"),]
    dfPI633731 <- B[which(B[,a]=="PI633731"),]
    dfPI639283 <- B[which(B[,a]=="PI639283"),]
    dfPI639285 <- B[which(B[,a]=="PI639285"),]
    dfPI639693 <- B[which(B[,a]=="PI639693"),]
    dfPI639740 <- B[which(B[,a]=="PI639740"),]
    dfPI643146 <- B[which(B[,a]=="PI643146"),]
    dfPI643395 <- B[which(B[,a]=="PI643395"),]
    dfPI660989 <- B[which(B[,a]=="PI660989"),]
    dfPI664025 <- B[which(B[,a]=="PI664025"),]
    dfS0613640 <- B[which(B[,a]=="S06.13640"),]
    dfTN053027 <- B[which(B[,a]=="TN05.3027"),]
    dfU03100612 <- B[which(B[,a]=="U03.100612"),]
    dfX4J10534 <- B[which(B[,a]=="X4J105.3.4"),]
    dfX5M20252 <- B[which(B[,a]=="X5M20.2.5.2"),]
    aeematrix[1,a] <- mean(dfIA3023$aeematrix)-mean(B$aeematrix)
    aeematrix[2,a] <- mean(dfCL0J09546$aeematrix)-mean(B$aeematrix)
    aeematrix[3,a] <- mean(dfHS63976$aeematrix)-mean(B$aeematrix)
    aeematrix[4,a] <- mean(dfLD015907$aeematrix)-mean(B$aeematrix)
    aeematrix[5,a] <- mean(dfLD024485$aeematrix)-mean(B$aeematrix)
    aeematrix[6,a] <- mean(dfLD029050$aeematrix)-mean(B$aeematrix)
    aeematrix[7,a] <- mean(dfLG032979$aeematrix)-mean(B$aeematrix)
    aeematrix[8,a] <- mean(dfLG033191$aeematrix)-mean(B$aeematrix)
    aeematrix[9,a] <- mean(dfLG044717$aeematrix)-mean(B$aeematrix)
    aeematrix[10,a] <- mean(dfLG054292$aeematrix)-mean(B$aeematrix)
    aeematrix[11,a] <- mean(dfLG054317$aeematrix)-mean(B$aeematrix)
    aeematrix[12,a] <- mean(dfLG054464$aeematrix)-mean(B$aeematrix)
    aeematrix[13,a] <- mean(dfLG054832$aeematrix)-mean(B$aeematrix)
    aeematrix[14,a] <- mean(dfNE3001$aeematrix)-mean(B$aeematrix)
    aeematrix[15,a] <- mean(dfPI398881$aeematrix)-mean(B$aeematrix)
    aeematrix[16,a] <- mean(dfPI404188A$aeematrix)-mean(B$aeematrix)
    aeematrix[17,a] <- mean(dfPI427136$aeematrix)-mean(B$aeematrix)
    aeematrix[18,a] <- mean(dfPI437169B$aeematrix)-mean(B$aeematrix)
    aeematrix[19,a] <- mean(dfPI518751$aeematrix)-mean(B$aeematrix)
    aeematrix[20,a] <- mean(dfPI561370$aeematrix)-mean(B$aeematrix)
    aeematrix[21,a] <- mean(dfPI574486$aeematrix)-mean(B$aeematrix)
    aeematrix[22,a] <- mean(dfPI595362$aeematrix)-mean(B$aeematrix)
    aeematrix[23,a] <- mean(dfPI598124$aeematrix)-mean(B$aeematrix)
    aeematrix[24,a] <- mean(dfPI602995$aeematrix)-mean(B$aeematrix)
    aeematrix[25,a] <- mean(dfPI615553$aeematrix)-mean(B$aeematrix)
    aeematrix[26,a] <- mean(dfPI633730$aeematrix)-mean(B$aeematrix)
    aeematrix[27,a] <- mean(dfPI633731$aeematrix)-mean(B$aeematrix)
    aeematrix[28,a] <- mean(dfPI639283$aeematrix)-mean(B$aeematrix)
    aeematrix[29,a] <- mean(dfPI639285$aeematrix)-mean(B$aeematrix)
    aeematrix[30,a] <- mean(dfPI639693$aeematrix)-mean(B$aeematrix)
    aeematrix[31,a] <- mean(dfPI639740$aeematrix)-mean(B$aeematrix)
    aeematrix[32,a] <- mean(dfPI643146$aeematrix)-mean(B$aeematrix)
    aeematrix[33,a] <- mean(dfPI643395$aeematrix)-mean(B$aeematrix)
    aeematrix[34,a] <- mean(dfPI660989$aeematrix)-mean(B$aeematrix)
    aeematrix[35,a] <- mean(dfPI664025$aeematrix)-mean(B$aeematrix)
    aeematrix[36,a] <- mean(dfS0613640$aeematrix)-mean(B$aeematrix)
    aeematrix[37,a] <- mean(dfTN053027$aeematrix)-mean(B$aeematrix)
    aeematrix[38,a] <- mean(dfU03100612$aeematrix)-mean(B$aeematrix)
    aeematrix[39,a] <- mean(dfX4J10534$aeematrix)-mean(B$aeematrix)
    aeematrix[40,a] <- mean(dfX5M20252$aeematrix)-mean(B$aeematrix)
    aeematrix[aeematrix=="NaN"] <- 0
    aeematrix[41,a] <- 0.5*(aeematrix[1,a]+aeematrix[2,a])
    aeematrix[42,a] <- 0.5*(aeematrix[1,a]+aeematrix[3,a])
    aeematrix[43,a] <- 0.5*(aeematrix[1,a]+aeematrix[4,a])
    aeematrix[44,a] <- 0.5*(aeematrix[1,a]+aeematrix[5,a])
    aeematrix[45,a] <- 0.5*(aeematrix[1,a]+aeematrix[6,a])
    aeematrix[46,a] <- 0.5*(aeematrix[1,a]+aeematrix[7,a])
    aeematrix[47,a] <- 0.5*(aeematrix[1,a]+aeematrix[8,a])
    aeematrix[48,a] <- 0.5*(aeematrix[1,a]+aeematrix[9,a])
    aeematrix[49,a] <- 0.5*(aeematrix[1,a]+aeematrix[10,a])
    aeematrix[50,a] <- 0.5*(aeematrix[1,a]+aeematrix[11,a])
    aeematrix[51,a] <- 0.5*(aeematrix[1,a]+aeematrix[12,a])
    aeematrix[52,a] <- 0.5*(aeematrix[1,a]+aeematrix[13,a])
    aeematrix[53,a] <- 0.5*(aeematrix[1,a]+aeematrix[14,a])
    aeematrix[54,a] <- 0.5*(aeematrix[1,a]+aeematrix[15,a])
    aeematrix[55,a] <- 0.5*(aeematrix[1,a]+aeematrix[16,a])
    aeematrix[56,a] <- 0.5*(aeematrix[1,a]+aeematrix[17,a])
    aeematrix[57,a] <- 0.5*(aeematrix[1,a]+aeematrix[18,a])
    aeematrix[58,a] <- 0.5*(aeematrix[1,a]+aeematrix[19,a])
    aeematrix[59,a] <- 0.5*(aeematrix[1,a]+aeematrix[20,a])
    aeematrix[60,a] <- 0.5*(aeematrix[1,a]+aeematrix[21,a])
    aeematrix[61,a] <- 0.5*(aeematrix[1,a]+aeematrix[22,a])
    aeematrix[62,a] <- 0.5*(aeematrix[1,a]+aeematrix[23,a])
    aeematrix[63,a] <- 0.5*(aeematrix[1,a]+aeematrix[24,a])
    aeematrix[64,a] <- 0.5*(aeematrix[1,a]+aeematrix[25,a])
    aeematrix[65,a] <- 0.5*(aeematrix[1,a]+aeematrix[26,a])
    aeematrix[66,a] <- 0.5*(aeematrix[1,a]+aeematrix[27,a])
    aeematrix[67,a] <- 0.5*(aeematrix[1,a]+aeematrix[28,a])
    aeematrix[68,a] <- 0.5*(aeematrix[1,a]+aeematrix[29,a])
    aeematrix[69,a] <- 0.5*(aeematrix[1,a]+aeematrix[30,a])
    aeematrix[70,a] <- 0.5*(aeematrix[1,a]+aeematrix[31,a])
    aeematrix[71,a] <- 0.5*(aeematrix[1,a]+aeematrix[32,a])
    aeematrix[72,a] <- 0.5*(aeematrix[1,a]+aeematrix[33,a])
    aeematrix[73,a] <- 0.5*(aeematrix[1,a]+aeematrix[34,a])
    aeematrix[74,a] <- 0.5*(aeematrix[1,a]+aeematrix[35,a])
    aeematrix[75,a] <- 0.5*(aeematrix[1,a]+aeematrix[36,a])
    aeematrix[76,a] <- 0.5*(aeematrix[1,a]+aeematrix[37,a])
    aeematrix[77,a] <- 0.5*(aeematrix[1,a]+aeematrix[38,a])
    aeematrix[78,a] <- 0.5*(aeematrix[1,a]+aeematrix[39,a])
    aeematrix[79,a] <- 0.5*(aeematrix[1,a]+aeematrix[40,a])
    replacewith <- as.data.frame(aeematrix[1:79,a])
    colnames(replacewith) <- "AEE"
    replacewith$Lookup <- c("IA3023","CL0J095.4.6","HS6.3976","LD01.5907","LD02.4485","LD02.9050","LG03.2979","LG03.3191","LG04.4717","LG05.4292","LG05.4317","LG05.4464","LG05.4832","NE3001","PI398881","PI404188A","PI427136","PI437169B","PI518751","PI561370","PI574486","PI595362","PI598124","PI602995","PI615553","PI633730","PI633731","PI639283","PI639285","PI639693","PI639740","PI643146","PI643395","PI660989","PI664025","S06.13640","TN05.3027","U03.100612","X4J105.3.4","X5M20.2.5.2","IA3023,CL0J095.4.6","HS6.3976,IA3023","IA3023,LD01.5907","IA3023,LD02.4485","IA3023,LD02.9050","IA3023,LG03.2979","IA3023,LG03.3191","IA3023,LG04.4717","IA3023,LG05.4292","IA3023,LG05.4317","IA3023,LG05.4464","IA3023,LG05.4832","IA3023,NE3001","IA3023,PI398881","IA3023,PI404188A","IA3023,PI427136","IA3023,PI437169B","IA3023,PI518751","IA3023,PI561370","IA3023,PI574486","IA3023,PI595362","IA3023,PI598124","IA3023,PI602995","IA3023,PI615553","IA3023,PI633730","IA3023,PI633731","IA3023,PI639283","IA3023,PI639285","IA3023,PI639693","IA3023,PI639740","IA3023,PI643146","IA3023,PI643395","IA3023,PI660989","IA3023,PI664025","IA3023,S06.13640","TN05.3027,IA3023","IA3023,U03.100612","IA3023,X4J105.3.4","IA3023,X5M20.2.5.2")
    replacement <- replacewith[seq(dim(replacewith)[1],1),]
    SNP_training_data[,a] <- stri_replace_all_fixed(SNP_training_data[,a], pattern = replacement$Lookup, replacement = replacement$AEE, vectorize_all = FALSE)
    SNP_testing_data[,a] <- stri_replace_all_fixed(SNP_testing_data[,a], pattern = replacement$Lookup, replacement = replacement$AEE, vectorize_all = FALSE)
  }
  SNP_testing_data <- data.frame(SNP_testing_data)
  SNP_training_data <- data.frame(SNP_training_data)
  fwrite(SNP_training_data,"train.temp")
  b1 <- fread("train.temp",colClasses = "numeric")
  fwrite(SNP_testing_data,"test.temp")
  b2 <- fread("test.temp",colClasses = "numeric")
  trained_model <- mixed.solve(y=Pheno_training_data,Z=b1)
  marker_effects <- as.matrix(trained_model$u)
  BLUE <- as.vector(trained_model$beta)
  predicted_train = as.matrix(b1) %*% marker_effects
  predicted_test = as.matrix(b2) %*% marker_effects
  predicted_train_result <- as.vector((predicted_train[,1])+BLUE)
  predicted_test_result <- as.vector((predicted_test[,1])+BLUE)
  trained_model2 <- mixed.solve(y=Pheno_training_data,Z=SNP_training_data2)
  marker_effects2 <- as.matrix(trained_model2$u)
  BLUE2 <- as.vector(trained_model2$beta)
  predicted_train2 = as.matrix(SNP_training_data2) %*% marker_effects2
  predicted_test2 = as.matrix(SNP_testing_data2) %*% marker_effects2
  predicted_train_result2 <- as.vector((predicted_train2[,1])+BLUE2)
  predicted_test_result2 <- as.vector((predicted_test2[,1])+BLUE2)
  Recordcorrelations[q,1] <- p    
  Recordcorrelations[q,2] <- cor(as.vector(Pheno_testing_data),predicted_test_result,use="complete")    
  Recordcorrelations[q,3] <- cor(as.vector(Pheno_training_data),predicted_train_result,use="complete")
  Recordcorrelations[q,4] <- cor(as.vector(Pheno_testing_data),predicted_test_result2,use="complete")    
  Recordcorrelations[q,5] <- cor(as.vector(Pheno_training_data),predicted_train_result2,use="complete")
}
colnames(Recordcorrelations) <- c("Environment","TestingAEE","TrainingAEE","TestingStandard","TrainingStandard")
view(Recordcorrelations)
write.csv(Recordcorrelations,'Standard vs AEEDTM.csv')
