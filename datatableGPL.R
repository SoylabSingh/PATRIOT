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
options(datatable.WhenJisSymbolThenCallingScope=TRUE)
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
  s <- s %>% drop_na(Yldkgha) #change based on trait
  df1 <- s[which(rownames(s) %in% rownames(A1)),]
  A2 <- A1[which(rownames(A1) %in% rownames(df1)),]
  df2 <- s[which(rownames(s) %in% rownames(A3)),]
  A4 <- A3[which(rownames(A3) %in% rownames(df2)),]
  u <- nrow(df1)
  Pheno <- df1 %>% select(Yldkgha) #change based on trait
  training_entries <- as.matrix(sample(1:u,floor(u*0.8)))
  testing_entries <- setdiff(1:u,training_entries)
  Pheno_training_data = as.matrix(Pheno[training_entries,])
  colnames(Pheno_training_data) <- "aeematrix"
  Pheno_testing_data = as.matrix(Pheno[testing_entries,])
  SNP_training_data = as.data.frame(A2[training_entries,])
  SNP_testing_data = as.data.frame(A2[testing_entries,])  
  SNP_training_data2 = as.data.frame(A4[training_entries,])
  SNP_testing_data2 = as.data.frame(A4[testing_entries,])
  B <- as.data.frame(cbind(SNP_training_data,Pheno_training_data))
  aeematrix <- data.frame()
  for (a in 1:(ncol(B)-1)) {
    B1 <- as.data.table(cbind(as.data.frame(B[,a]),B[,4290]))
    colnames(B1) <- c("a","b")
    aeematrix[1,a] <- B1["IA3023",mean(b),on="a"]-mean(B1$b)
    aeematrix[2,a] <- B1["CL0J095.4.6",mean(b),on="a"]-mean(B1$b)
    aeematrix[3,a] <- B1["HS6.3976",mean(b),on="a"]-mean(B1$b)
    aeematrix[4,a] <- B1["LD01.5907",mean(b),on="a"]-mean(B1$b)
    aeematrix[5,a] <- B1["LD02.4485",mean(b),on="a"]-mean(B1$b)
    aeematrix[6,a] <- B1["LD02.9050",mean(b),on="a"]-mean(B1$b)
    aeematrix[7,a] <- B1["LG03.2979",mean(b),on="a"]-mean(B1$b)
    aeematrix[8,a] <- B1["LG03.3191",mean(b),on="a"]-mean(B1$b)
    aeematrix[9,a] <- B1["LG04.4717",mean(b),on="a"]-mean(B1$b)
    aeematrix[10,a] <- B1["LG05.4292",mean(b),on="a"]-mean(B1$b)
    aeematrix[11,a] <- B1["LG05.4317",mean(b),on="a"]-mean(B1$b)
    aeematrix[12,a] <- B1["LG05.4464",mean(b),on="a"]-mean(B1$b)
    aeematrix[13,a] <- B1["LG05.4832",mean(b),on="a"]-mean(B1$b)
    aeematrix[14,a] <- B1["NE3001",mean(b),on="a"]-mean(B1$b)
    aeematrix[15,a] <- B1["PI398881",mean(b),on="a"]-mean(B1$b)
    aeematrix[16,a] <- B1["PI404188A",mean(b),on="a"]-mean(B1$b)
    aeematrix[17,a] <- B1["PI427136",mean(b),on="a"]-mean(B1$b)
    aeematrix[18,a] <- B1["PI437169B",mean(b),on="a"]-mean(B1$b)
    aeematrix[19,a] <- B1["PI518751",mean(b),on="a"]-mean(B1$b)
    aeematrix[20,a] <- B1["PI561370",mean(b),on="a"]-mean(B1$b)
    aeematrix[21,a] <- B1["PI574486",mean(b),on="a"]-mean(B1$b)
    aeematrix[22,a] <- B1["PI595362",mean(b),on="a"]-mean(B1$b)
    aeematrix[23,a] <- B1["PI598124",mean(b),on="a"]-mean(B1$b)
    aeematrix[24,a] <- B1["PI602995",mean(b),on="a"]-mean(B1$b)
    aeematrix[25,a] <- B1["PI615553",mean(b),on="a"]-mean(B1$b)
    aeematrix[26,a] <- B1["PI633730",mean(b),on="a"]-mean(B1$b)
    aeematrix[27,a] <- B1["PI633731",mean(b),on="a"]-mean(B1$b)
    aeematrix[28,a] <- B1["PI639283",mean(b),on="a"]-mean(B1$b)
    aeematrix[29,a] <- B1["PI639285",mean(b),on="a"]-mean(B1$b)
    aeematrix[30,a] <- B1["PI639693",mean(b),on="a"]-mean(B1$b)
    aeematrix[31,a] <- B1["PI639740",mean(b),on="a"]-mean(B1$b)
    aeematrix[32,a] <- B1["PI643146",mean(b),on="a"]-mean(B1$b)
    aeematrix[33,a] <- B1["PI643395",mean(b),on="a"]-mean(B1$b)
    aeematrix[34,a] <- B1["PI660989",mean(b),on="a"]-mean(B1$b)
    aeematrix[35,a] <- B1["PI664025",mean(b),on="a"]-mean(B1$b)
    aeematrix[36,a] <- B1["S06.13640",mean(b),on="a"]-mean(B1$b)
    aeematrix[37,a] <- B1["TN05.3027",mean(b),on="a"]-mean(B1$b)
    aeematrix[38,a] <- B1["U03.100612",mean(b),on="a"]-mean(B1$b)
    aeematrix[39,a] <- B1["X4J105.3.4",mean(b),on="a"]-mean(B1$b)
    aeematrix[40,a] <- B1["X5M20.2.5.2",mean(b),on="a"]-mean(B1$b)
    aeematrix[is.na(aeematrix)] <- 0
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
write.csv(Recordcorrelations,'Standard vs AEEyld.csv')