#!/usr/bin/env Rscript

### Script to run resampling stability analysis

library(corrplot)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape)

insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

subjects <- read.csv("data/Data_for_SNF/ids.csv", header=FALSE)
CT <- read.csv("data/Data_for_SNF/CT_snf.csv", header=FALSE)
volume <- read.csv("data/Data_for_SNF/Volumes_snf.csv", header=FALSE)
FA <- read.csv("data/Data_for_SNF/FA_snf.csv", header=FALSE)
clinical <-  read.csv("data/Data_for_SNF/pond_snf.csv", header=FALSE)

directory <- ("output/3_Stability")

source("code/3_Stability_testing/bootstrapping_scripts/bootstrapping.r")
numboot=1000
nsub=176
K=18
alpha=0.8
t=10
bootsize=0.8
clusters=4

# Setting up which participants will be included in each permutation
permutation_matrix <- bootstrapping_SNF(numboot=numboot, nsub=nsub, bootsize=bootsize)
# Getting clustering solutions for all the permuatations of sampled participants using SNF 
clus_sil <- clustering(perms=permutation_matrix, bootsize=bootsize, K=K, t=t, alpha=alpha, clusters=clusters, CT=CT, FA=FA, volume=volume, clinical=clinical)
# Dividing output matrix into clusters and silhouette widths
clus_out <- clus_sil[1:(numboot), ]
silhouette_width <- clus_sil[(numboot+1):(numboot*2), ]
## getting NMI scores for each permutation
All_NMI_scores <- NMI_scores(perms=permutation_matrix, bootsize=bootsize, K=K, t=t, alpha=alpha, clusters=clusters, CT=CT, FA=FA, clinical=clinical, volume=volume)
# getting the adjusted rand index between all clustering solutions
list_randindex <- stability(clus_out=clus_out, perms=permutation_matrix) # returns brandindex:adjusted rand index
# Calculate how often each participant is clusted together and the probability that they will be clustered together
percent_agree <- percent_agree(clus_out=clus_out)

#list_adjustedrandindex <- list_randindex[ ,length(numboot):1000]
# list_randindex_new <- list_randindex[,(1:500)]

write.csv(percent_agree, file=file.path(directory, paste("Percent_agree_4c_1000perms.csv", sep="")))
write.csv(list_randindex, file=file.path(directory, paste("Rand_indices_4c_1000perms.csv", sep="")))

write.csv(clus_out, file=file.path(directory, paste("Adj_rand_indices_4c_1000perms.csv", sep="")))
write.csv(permutation_matrix, file=file.path(directory, paste("permutation_matrix_4c_1000perms.csv", sep="")))

write.csv(All_NMI_scores, file=file.path(directory, paste("NMI_scores_0.8_1000perms.csv", sep="")))

##### Sorting out measures and getting the top 35
measures <- read.csv("data/measure_names.csv")
NMI <- read.csv(file=file.path(directory, paste("NMI_scores_0.8_1000perms.csv", sep="")))
NMI$X <- NULL

top_NMI_measures <- data.frame(matrix(0, nrow = 35, ncol = numboot))
i=1
for (i in seq(1:length(NMI))) {
  intermediate <- NMI[ , i]
  intermediate <- cbind(measures, intermediate)
  intermediate$intermediate <- as.numeric(as.character(intermediate$intermediate))
  intermediate <- intermediate[order(-intermediate$intermediate), ]

  top_NMI_measures[ ,i] <- intermediate[1:35, 1]
}

write.csv(top_NMI_measures, file=file.path(directory, paste("Top_35_NMI_scores_0.8_1000perms.csv", sep="")))

#### counting how often measures are in the top 35
top_measures <- top_NMI_measures
top_measures$Column <- seq(1:length(top_measures[,1]))
top_measures <- melt(top_measures, id.vars=c('Column'),var='Index')
top_measures$Column <- NULL
top_measures$Index <- NULL

measure_count <- as.data.frame(table(top_measures$value))
names(measure_count)[names(measure_count) == "Var1"] <- "Measure"
measure_count <- measure_count[order(measure_count$Freq),]

# add in that there are 0 times IFO-L is in the top 35 measures
write.csv(measure_count, file=file.path(directory, paste("Count-top_35_NMI_0.8_1000perms.csv", sep="")))

