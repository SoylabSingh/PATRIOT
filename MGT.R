setwd('/Users/User/Desktop/Imputed/')
df <- read.csv('868ti.csv')
df <- df[1:42080,1:871]
df2 <- as.matrix(df)

for (x in 1:42080){
  for (y in 4:871){
    c <- as.character(df[x,y])
    if(c %in% colnames(df)){
      beep <- df[c]
      df2[x,y] <- as.character(beep[x,1])
    }
  }
}
write.csv(df2,'MGT2.csv')