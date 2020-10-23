install.packages("snow")
install.packages("doSNOW")
install.packages("parallel")
install.packages("rrBLUP")
install.packages("dplyr")
library(snow)
library(parallel)
library(doSNOW)
library(ggplot2)
detectCores()
cl <- makeCluster(16,type="SOCK")
registerDoSNOW(cl)
library(rrBLUP)
library(dplyr)
setwd('/Users/User/Desktop/Imputed')
Pheno_data <- read.csv('NAMFullPheno.csv')

Pheno <- data_2 %>% select(Yldkgha)

GWAS_GD <- read.csv("NAMt.csv")
rownames(GWAS_GD) <- GWAS_GD$X
GWAS_GD$X <- NULL
GWAS_GD[1:5,1:5]
data_2 <- read.csv("NAMPhenoLD0033092013Diers.csv")
GD2 <- read.csv("NAMLD003309.csv")
rownames(GD2) <- GD2$X
GD2$X <- NULL
GD2 <- GD2[sort(rownames(GD2)),]
GD1 <- GWAS_GD[sort(rownames(GWAS_GD)),]
GD1 <- GD1[which(rownames(GD1) %in% data_2$CorrectedStrain),]

training_entries <- as.matrix(sample(1:141,100))
testing_entries <- setdiff(1:141,training_entries)
Pheno_training_data = as.matrix(Pheno[training_entries,])
Pheno_testing_data = as.matrix(Pheno[testing_entries,])
SNP_training_data = as.matrix(GD1[training_entries,])
SNP_testing_data = as.matrix(GD1[testing_entries,])

trained_model <- mixed.solve(y=Pheno_training_data,Z=SNP_training_data)
marker_effects <- as.matrix(trained_model$u)
BLUE <- as.vector(trained_model$beta)
predicted_train = as.matrix(SNP_training_data) %*% marker_effects
predicted_test = as.matrix(SNP_testing_data) %*% marker_effects
predicted_train_result <- as.vector((predicted_train[,1])+BLUE)
predicted_test_result <- as.vector((predicted_test[,1])+BLUE)
cor(as.vector(Pheno_testing_data),predicted_test_result,use="complete")
cor(as.vector(Pheno_training_data),predicted_train_result,use="complete")

reg <- (lm(predicted_train_result ~ Pheno_training_data))
coeff=coefficients(reg)
eq = paste0("y = ", round(coeff[2],1), "x+", round(coeff[1],1))
plot(Pheno_training_data,predicted_train_result,main=eq,ylab="Prediction without tracing")
abline(lm(predicted_train_result ~ Pheno_training_data))


training_entries <- as.matrix(sample(1:141,100))
testing_entries <- setdiff(1:141,training_entries)
SNP_training_data2 = as.matrix(GD2[training_entries,],K=NULL)
SNP_testing_data2 = as.matrix(GD2[testing_entries,],K=NULL)
trained_model2 <- mixed.solve(y=Pheno_training_data,Z=SNP_training_data2)
marker_effects2 <- as.matrix(trained_model2$u)
BLUE2 <- as.vector(trained_model2$beta)
predicted_train2 = as.matrix(SNP_training_data2) %*% marker_effects2
predicted_test2 = as.matrix(SNP_testing_data2) %*% marker_effects2
predicted_train_result2 <- as.vector((predicted_train2[,1])+BLUE2)
predicted_test_result2 <- as.vector((predicted_test2[,1])+BLUE2)
cor(as.vector(Pheno_testing_data),predicted_test_result2,use="complete")
cor(as.vector(Pheno_training_data),predicted_train_result2,use="complete")
reg2 <- (lm(predicted_train_result2 ~ Pheno_training_data))
coeff2=coefficients(reg2)
eq2 = paste0("y = ", round(coeff2[2],1), "x+", round(coeff2[1],1))
plot(Pheno_training_data,predicted_train_result2,main=eq2,ylab="Prediction with tracing")
abline(lm(predicted_train_result2 ~ Pheno_training_data))
#####Using raw marker data#####
df <- Pheno_data
df <- df[which(df$Env=="2012_IA"),]
markers <- read.csv('NAMt.csv')
rownames(markers) <- markers$X
markers$X <- NULL
markers$CHROM <- NULL
markers$POS <- NULL
df <- df1[which(df1$CorrectedStrain %in% rownames(markers)),]
markers1 <- markers[which(rownames(markers) %in% df$CorrectedStrain),]
Pheno <- df %>% select(Yldkgha)

training_entries <- as.matrix(sample(1:5111,4000))
testing_entries <- setdiff(1:5111,training_entries)
Pheno_training_data = as.matrix(Pheno[training_entries,])
Pheno_testing_data = as.matrix(Pheno[testing_entries,])
SNP_training_data = as.matrix(markers1[training_entries,])
SNP_testing_data = as.matrix(markers1[testing_entries,])

trained_model <- mixed.solve(y=Pheno_training_data[,11],Z=beta)
marker_effects <- as.matrix(trained_model$u)
BLUE <- as.vector(trained_model$beta)
predicted_train = as.matrix(SNP_training_data) %*% marker_effects
predicted_test = as.matrix(SNP_testing_data) %*% marker_effects
predicted_train_result <- as.vector((predicted_train[,1])+BLUE)
predicted_test_result <- as.vector((predicted_test[,1])+BLUE)
cor(as.vector(Pheno_testing_data),predicted_test_result,use="complete")
cor(as.vector(Pheno_training_data),predicted_train_result,use="complete")
######Using Chromosome Tracing and allele effect estimates#####
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
df1 <- df[which(df$Env=="2012_IA"),]
rownames(df1) <- df1$CorrectedStrain
df <- df1[which(df1$CorrectedStrain %in% rownames(A1)),]
df2 <- df[which(df$Yldkgha> 0),]
A <- A1[which(rownames(A1) %in% df2$CorrectedStrain),]
A2 <- A3[which(rownames(A3) %in% df2$CorrectedStrain),]
u <- ncol(A2)
training_entries <- as.matrix(sample(1:u,0.8*u))
testing_entries <- setdiff(1:u,training_entries)
Pheno_training_data = as.matrix(df2[training_entries,])
Pheno_testing_data = as.matrix(df2[testing_entries,])
SNP_training_data = as.matrix(A[training_entries,])
SNP_testing_data = as.matrix(A[testing_entries,])
SNP_training_data1 = as.matrix(A2[training_entries,])
SNP_testing_data1 = as.matrix(A2[testing_entries,])

#####
B <- cbind(A,df2)
Yldkgha <- data.frame()
for (a in 1:(ncol(B)-16)) {
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
  Yldkgha[1,a] <- mean(dfIA3023$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[2,a] <- mean(dfCL0J09546$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[3,a] <- mean(dfHS63976$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[4,a] <- mean(dfLD015907$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[5,a] <- mean(dfLD024485$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[6,a] <- mean(dfLD029050$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[7,a] <- mean(dfLG032979$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[8,a] <- mean(dfLG033191$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[9,a] <- mean(dfLG044717$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[10,a] <- mean(dfLG054292$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[11,a] <- mean(dfLG054317$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[12,a] <- mean(dfLG054464$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[13,a] <- mean(dfLG054832$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[14,a] <- mean(dfNE3001$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[15,a] <- mean(dfPI398881$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[16,a] <- mean(dfPI404188A$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[17,a] <- mean(dfPI427136$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[18,a] <- mean(dfPI437169B$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[19,a] <- mean(dfPI518751$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[20,a] <- mean(dfPI561370$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[21,a] <- mean(dfPI574486$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[22,a] <- mean(dfPI595362$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[23,a] <- mean(dfPI598124$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[24,a] <- mean(dfPI602995$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[25,a] <- mean(dfPI615553$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[26,a] <- mean(dfPI633730$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[27,a] <- mean(dfPI633731$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[28,a] <- mean(dfPI639283$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[29,a] <- mean(dfPI639285$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[30,a] <- mean(dfPI639693$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[31,a] <- mean(dfPI639740$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[32,a] <- mean(dfPI643146$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[33,a] <- mean(dfPI643395$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[34,a] <- mean(dfPI660989$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[35,a] <- mean(dfPI664025$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[36,a] <- mean(dfS0613640$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[37,a] <- mean(dfTN053027$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[38,a] <- mean(dfU03100612$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[39,a] <- mean(dfX4J10534$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[40,a] <- mean(dfX5M20252$Yldkgha)-mean(B$Yldkgha)
  Yldkgha[Yldkgha=="NaN"] <- 0
  Yldkgha[41,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[2,a])
  Yldkgha[42,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[3,a])
  Yldkgha[43,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[4,a])
  Yldkgha[44,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[5,a])
  Yldkgha[45,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[6,a])
  Yldkgha[46,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[7,a])
  Yldkgha[47,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[8,a])
  Yldkgha[48,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[9,a])
  Yldkgha[49,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[10,a])
  Yldkgha[50,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[11,a])
  Yldkgha[51,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[12,a])
  Yldkgha[52,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[13,a])
  Yldkgha[53,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[14,a])
  Yldkgha[54,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[15,a])
  Yldkgha[55,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[16,a])
  Yldkgha[56,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[17,a])
  Yldkgha[57,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[18,a])
  Yldkgha[58,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[19,a])
  Yldkgha[59,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[20,a])
  Yldkgha[60,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[21,a])
  Yldkgha[61,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[22,a])
  Yldkgha[62,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[23,a])
  Yldkgha[63,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[24,a])
  Yldkgha[64,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[25,a])
  Yldkgha[65,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[26,a])
  Yldkgha[66,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[27,a])
  Yldkgha[67,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[28,a])
  Yldkgha[68,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[29,a])
  Yldkgha[69,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[30,a])
  Yldkgha[70,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[31,a])
  Yldkgha[71,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[32,a])
  Yldkgha[72,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[33,a])
  Yldkgha[73,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[34,a])
  Yldkgha[74,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[35,a])
  Yldkgha[75,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[36,a])
  Yldkgha[76,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[37,a])
  Yldkgha[77,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[38,a])
  Yldkgha[78,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[39,a])
  Yldkgha[79,a] <- 0.5*(Yldkgha[1,a]+Yldkgha[40,a])
  replacewith <- as.data.frame(Yldkgha[1:79,a])
  colnames(replacewith) <- "AEE"
  replacewith$Lookup <- c("IA3023","CL0J095.4.6","HS6.3976","LD01.5907","LD02.4485","LD02.9050","LG03.2979","LG03.3191","LG04.4717","LG05.4292","LG05.4317","LG05.4464","LG05.4832","NE3001","PI398881","PI404188A","PI427136","PI437169B","PI518751","PI561370","PI574486","PI595362","PI598124","PI602995","PI615553","PI633730","PI633731","PI639283","PI639285","PI639693","PI639740","PI643146","PI643395","PI660989","PI664025","S06.13640","TN05.3027","U03.100612","X4J105.3.4","X5M20.2.5.2","IA3023,CL0J095.4.6","HS6.3976,IA3023","IA3023,LD01.5907","IA3023,LD02.4485","IA3023,LD02.9050","IA3023,LG03.2979","IA3023,LG03.3191","IA3023,LG04.4717","IA3023,LG05.4292","IA3023,LG05.4317","IA3023,LG05.4464","IA3023,LG05.4832","IA3023,NE3001","IA3023,PI398881","IA3023,PI404188A","IA3023,PI427136","IA3023,PI437169B","IA3023,PI518751","IA3023,PI561370","IA3023,PI574486","IA3023,PI595362","IA3023,PI598124","IA3023,PI602995","IA3023,PI615553","IA3023,PI633730","IA3023,PI633731","IA3023,PI639283","IA3023,PI639285","IA3023,PI639693","IA3023,PI639740","IA3023,PI643146","IA3023,PI643395","IA3023,PI660989","IA3023,PI664025","IA3023,S06.13640","TN05.3027,IA3023","IA3023,U03.100612","IA3023,X4J105.3.4","IA3023,X5M20.2.5.2")
  replacement <- replacewith[seq(dim(replacewith)[1],1),]
  SNP_training_data[,a] <- stri_replace_all_fixed(SNP_training_data[,a], pattern = replacement$Lookup, replacement = replacement$AEE, vectorize_all = FALSE)
  SNP_testing_data[,a] <- stri_replace_all_fixed(SNP_testing_data[,a], pattern = replacement$Lookup, replacement = replacement$AEE, vectorize_all = FALSE)
}
dim(SNP_training_data)
beta <- as.matrix(SNP_training_data)
class(beta) <- "numeric"
beta[1:5,1:5]
class(Pheno_training_data) <- "numeric"
