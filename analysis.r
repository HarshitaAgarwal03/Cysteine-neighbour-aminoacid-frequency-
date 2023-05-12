#input file preparation commands in bash 
# $cut -f 1 mass-spec.tsv | tail -n 586 > protein_seq.csv
# $sed 's/.*].//g' protein_seq.csv | sed 's/.\[.*//g' > protein_string.csv
# $sed 's/C/_/g' protein_string.csv > protein_seq_input.csv



library(tidyverse)
library(ggplot2)
library(stringr)

protein = read.csv("protein_seq_single_occurance_C.csv",header = FALSE)       #read the file into a dataframe
#head(protein)

protein_split = data.frame(do.call('rbind', strsplit(as.character(protein$V1),'_',fixed=TRUE))) #split the string based on the delimeter, which was added using bash commands
#head(protein_split)

protein_split$before = str_sub(protein_split$X1, start=-3)    #get 3 amino acid before the desired amino acid
#head(protein_split)

protein_split$after = str_sub(protein_split$X2, end=3)        #get 3 amino acid after the desired amino acid
#head(protein_split, 20)

#----------------------------------------------------------------------------------------------------------------
minus_string = strsplit(protein_split$before,split="")
#head(minus_string,20)
#length(minus_string)

for (m in 1:509)
  {
    n = minus_string[m][1]
    if ((length(unlist(n))) == 2)
      {
        da = c("X" ,(unlist(n)))
        minus_string[m][1] = c(list(da))
      }
  
    if ((length(unlist(n))) == 1)
      {
        db = c("X" , "X" ,(unlist(n)))
        minus_string[m][1] = c(list(db))
      }
  
    if ((length(unlist(n))) == 0)
      {
        dd = c("X" , "X" , "X")
        minus_string[m][1] = c(list(dd))
      }
  }

data = as.data.frame(do.call(rbind, minus_string))
colnames(data) = c("minus3" , "minus2" , "minus1")
#head(data,20)
#head(minus_string)

#--------------------------------------------------------------------------------------------------------------------------------------

plus_string = strsplit(protein_split$after,split="")
head(plus_string,20)
#length(minus_string)

for (o in 1:509)
{
  p = plus_string[o][1]
  if ((length(unlist(p))) == 2)
  {
    de = c((unlist(p)), "X")
    plus_string[o][1] = c(list(de))
  }
  
  if ((length(unlist(p))) == 1)
  {
    df = c((unlist(p)), "X" , "X")
    plus_string[o][1] = c(list(df))
  }
  
  if ((length(unlist(p))) == 0)
  {
    dg = c("X" , "X" , "X")
    plus_string[o][1] = c(list(dg))
  }
}

data1 = as.data.frame(do.call(rbind, plus_string))
colnames(data1) = c("plus1" , "plus2" , "plus3")
#head(data1,40)
#head(plus_string)

#-----------------------------------------------------------------------------------------------------------

protein_split = cbind(protein_split, data,data1)
#head(protein_split,40)

#------------------------------------------------------------------------------------------------------

minus3 = table(protein_split$minus3)
minus2 = table(protein_split$minus2)
minus1 = table(protein_split$minus1)
plus1 = table(protein_split$plus1)
plus2 = table(protein_split$plus2)
plus3 = table(protein_split$plus3)

count_matrix <- as.data.frame(rbind(minus3,minus2,minus1,plus1,plus2,plus3))

#-----------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
