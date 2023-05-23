#Code to calculate and plot frequency of amino acids present 3 position upstream and downstream of Cysteine and Dicysteine in a list of peptide sequences
#Classify neighbours of Cysteine and dicysteine as **Polar, non-polar, negatively charged and positively charged** in a list of peptide sequence

## Author: Harshita

# one could modify code according to individual's need

######--------------************************************************************************--------------######

setwd("~/Harshita/mass-spec-data")      #set working directory

#--------------------------------- loading packages for the analysis, if not present download using install.packages()------------------------

#install.packages("tidyverse")     #commands to install packages if not present
#install.packages("ggplot2")
#install.packages("stringr")
#install.packages("textworks")
#install.packages("dplyr")
#install.packages("RColorBrewer")

library(tidyverse)
library(ggplot2)
library(stringr)
library(textworks)
library(dplyr)
library(RColorBrewer)

#-------------------------------------------read data as data frame from file and split file into different types--------------------
#-------------------------------------------Bash commands run on the terminal-------------------------------------

# input file preparation commands in bash 
# input file name is protein_seq.csv
# $sed 's/.*].//g' protein_seq.csv | sed 's/.\[.*//g' > protein_string.csv
# $sed 's/C/_/g' protein_string.csv > protein_seq_input.csv
# $grep "__" protein_seq_input.csv > protein_seq_doubleC.csv
# $grep -v -E '_.*_|__' protein_seq_input.csv > protein_seq_single_occurance_C.csv
# $grep '_.*_' protein_seq_input.csv | grep -v "__" > protein_seq_multiple_occurance_C.csv

#--------------------------------------------------Algorithm-----------------------------------------------

protein = read.csv("protein_seq_single_occurance_C.csv",header = FALSE)       #read the file into a data frame 
#head(protein)                                                                #visualize dataframe
proteinD = read.csv("protein_seq_doubleC.csv",header = FALSE)       #read the file into a data frame
#head(proteinD)
proteinM = read.csv("protein_seq_multiple_occurance_C.csv",header = FALSE)  # read the file into a data frame 
#head(proteinM)
proteinm = rbind(proteinD,proteinM)   #combine data frame containing presence of >1 C

Nostrings = as.numeric(count(proteinm))   #count number of strings(protein sequences) with presence of >1 C 
#head(Nostrings)
amino = proteinm[['V1']]                  # save strings in a list
#head(amino)
amino_final = c()                         #declaring variables for downstream analysis 
a = list()
b = list()
c = list()
e = list()
f = list()
x = c()
j = 1
k = 1

#-------------*******************************************************************************************----------------------
#----------------------------------------------preparing protein sequences with >1 C for analysis-----------------------------
#--------------*******************************************************************************************----------------------

for (i in 1:Nostrings)
{
  x[i] = as.character(amino[i])
  if ((str_count(x[i],"_")) == 3)
  {
    a[[j]] = str_replace_nth(x = x[i],"_","C",n=3)
    a[[j]] = str_replace_nth(x = a[[j]], "_", "C",n=2)
    b[[j]] = str_replace_nth(x = x[i],"_","C",n=2)
    b[[j]] = str_replace_nth(x = b[[j]], "_", "C",n=1)
    c[[j]] = str_replace_nth(x = x[i],"_","C",n=3)
    c[[j]] = str_replace_nth(x = c[[j]], "_", "C", n=1)
    j = j+1
  }
  
  if ((str_count(x[i],"_")) == 2)
  {
    e[[k]] = str_replace_nth(x = x[i],"_","C",n=1)
    f[[k]] = str_replace_nth(x= x[i],"_","C",n=2)
    k = k+1
  }
}  

amino_final = c(unlist(a), unlist(b), unlist(c), unlist(e), unlist(f))           #create a list with updated sequences
amino_final_data = as.data.frame(amino_final)                                    # create a dataframe with the updated list
colnames(amino_final_data) = c("V1")                                             # assign column name to dataframe
protein = rbind(protein,amino_final_data)                                        #combine all the sequences

#-------------*****************************************************************************************************************-----------------------
#-----------------------------------make column with 3 amino acid before and after C in the sequences and save it in dataframe-----------------------
#-------------*****************************************************************************************************************----------------------

protein_split = data.frame(do.call('rbind', strsplit(as.character(protein$V1),'_',fixed=TRUE))) #split the string based on the delimeter, which was added using bash commands

protein_split$before = str_sub(protein_split$X1, start=-3)    #get 3 amino acid before the desired amino acid

protein_split$after = str_sub(protein_split$X2, end=3)        #get 3 amino acid after the desired amino acid

#----------********************************************************************************************************************-------------------------
#--------------------------------------add X to the positions where no amino acid is present----------------------------------------------------------
#----------********************************************************************************************************************-------------------------

minus_string = strsplit(protein_split$before,split="")


Nostrings_final = as.numeric(count(protein))

for (m in 1:Nostrings_final)
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

#-------------- positive positions------------------------------------

plus_string = strsplit(protein_split$after,split="")
head(plus_string,20)

for (o in 1:Nostrings_final)
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

#-----------------------------------------------------------------------------------------------------------

protein_split = cbind(protein_split, data,data1)

#-----------------------------Calculate no of amino acid present at different positions-----------------------

minus3 = table(protein_split$minus3)
minus2 = table(protein_split$minus2)
minus1 = table(protein_split$minus1)
plus1 = table(protein_split$plus1)
plus2 = table(protein_split$plus2)
plus3 = table(protein_split$plus3)

count_matrix <- as.data.frame(bind_rows(minus3,minus2,minus1,plus1,plus2,plus3))
row.names(count_matrix) <- c("minus3", "minus2","minus1","plus1","plus2","plus3")
head(count_matrix)
write.csv(count_matrix, file = "count_amino_acid.csv")         #save the results in a file

#---------------------------group amino acids into the group and count each groups number at each position---------------------------

Non_polar = c("G","A","V","C","P","L","I","M","W","F")         # group amino acids
polar = c("S","T","Y","N","Q")
plus_charge = c("K","R","H")
minus_charge = c("D","E")

non_polar_sum = rowSums(count_matrix[,Non_polar])
polar_sum = rowSums(count_matrix[,polar])
plus_charge_sum = rowSums(count_matrix[,plus_charge])
minus_charge_sum = rowSums(count_matrix[,minus_charge])

classify_aa = as.data.frame(bind_rows(non_polar_sum,polar_sum,plus_charge_sum,minus_charge_sum))
row.names(classify_aa) <- c("non_polar_AA","polar_AA","positive_charge_AA","negative_charge_AA")
head(classify_aa)
write.csv(classify_aa, file ="count_TypeOf_AA.csv")       #save the results in a file

#------------------------**********************************************************----------------
#-------------------------------------plot the counts-----------------------------------------------
#------------------------**********************************************************-------------------
count_matrix <- t(count_matrix)      #transpose count matrix for ploting barplots
count_matrix <- count_matrix[-20,]   #remove X empty coloum from the dataframe 
head(count_matrix)
graphcolor <- c("seashell","wheat","wheat2","pink","pink1","pink2","rosybrown1","rosybrown2","pink3","palevioletred","palevioletred1","palevioletred2","palevioletred2","palevioletred3","violetred1","violetred2","violetred3","violetred","palevioletred4","violetred4") #list of colors
jpeg("count_plot.jpg", width = 600, height = 350)  #save plot as jpg file
barplot(count_matrix, beside= T ,xlab="Position",ylab="Count of Amino acid ",width=0.2, col= graphcolor, legend.text = T, args.legend = list(x="topright",inset=c(-0.06,-0.1),bty="n"),main="Barplot of counts of Amino acid present upstream and downstream of Cysteine") #plot a barplot
dev.off()

classify_aa1 <- data.matrix(classify_aa)  #convert typesOfAA dataframe to matrix
coul <- brewer.pal(4, "Set2")             #set colours
jpeg("count_plotof_aa.jpg", width = 850, height = 750)  #save plot as jpg file
barplot(classify_aa1, beside= T,xlab="Position",ylab="Count of Amino acid",width=0.1, col= coul, legend.text = T, args.legend = list(x="topright",inset=c(-0.02,-0.03)),main="Barplot of counts of type of Amino acid present upstream and downstream of Cysteine") #plot a barplot
dev.off()

