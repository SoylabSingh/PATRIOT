#Set wd.#
setwd('/Users/User/Desktop/Imputed/')
#Read in marker file following tracing and imputation.#
df <- read.csv('868ti.csv')
#Subset df.#
df <- df[1:42080,1:871]
#Save df as matrix.#
df2 <- as.matrix(df)
#Loop through markers.#
for (x in 1:42080){
  #Loop through progeny. First three columns are skipped for Chrom, Pos, MarkerName.#
  for (y in 4:871){
    #Extract cell value.#
    c <- as.character(df[x,y])
    #If cell value is same as a column name...#
    if(c %in% colnames(df)){
      #Create df of that parent's traced calls.#
      beep <- df[c]
      #Overwrite initial cell value with corresponding value from the parent it inherited that SNP from.#
      df2[x,y] <- as.character(beep[x,1])}}}
#Write to file.#
write.csv(df2,'MGT2.csv')