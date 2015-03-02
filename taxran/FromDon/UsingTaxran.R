#this script shows how to use the function taxran() to calcualte the % variance
#attributible to each taxonomic level and simulte the null hypothesis of 
#no taxonomic structure to test the hypothesis of taxonomic conservatism.

#load data
setwd("/Users/don/Documents/DLM_Lennon")
dat=read.csv("nestedfinal.csv")
#load function
source("taxran15.R")

#peek at data
head(dat)

#run taxran for trait "Optimum"
ans=taxran15(dat,"Optimum")

#the returned object is a list with 2 slots:
#the pct var attributible to each tax level
ans[[1]]

#and a matrix of the results of the simulations
ans[[2]]

#find 2.5% and 97.5% percentile (i.e. 95% CI) of %attributible to Type given the null hypothesis
#actual
ans[[1]]["Type"]
#null CI
sort(ans[[2]]['Type',])[c(25,975)]
#conclusion: no reason to believe that there is sig conservation at the Type level

#How about at family level?
ans[[1]]["Family"]
sort(ans[[2]]['Family',])[c(25,975)]

#Wow. actual pct var attibutiable to family level is much greater than predicted by the 
#null hypothesis.
#conclusion: There is tax conservation at family level

#calc 1-sided pvalue
max(which(sort(ans[[2]]['Family',])>ans[[1]]['Family']),0)/1000
