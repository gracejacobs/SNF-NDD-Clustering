### Checking iterations of parameters and their effect on cluster number and NMI scores
library(dplyr)
library(SNFtool)

# creating output folders
dir.create("output/1_SNF_analysis/parameter_iterations/alpha_0.3")
dir.create("output/1_SNF_analysis/parameter_iterations/alpha_0.4")
dir.create("output/1_SNF_analysis/parameter_iterations/alpha_0.5")
dir.create("output/1_SNF_analysis/parameter_iterations/alpha_0.6")
dir.create("output/1_SNF_analysis/parameter_iterations/alpha_0.7")
dir.create("output/1_SNF_analysis/parameter_iterations/alpha_0.8")

# setting up empty dataframes to compare nearest neighbors (K) and alpha parameters
ktests=seq(10,30)
resultdf_1 = data.frame("K"=ktests,"0.3"=numeric(length(ktests)), "0.4"=numeric(length(ktests)), "0.5"=numeric(length(ktests)), "0.6"=numeric(length(ktests)), "0.7"=numeric(length(ktests)), "0.8"=numeric(length(ktests)))
resultdf_2 = data.frame("K"=ktests,"0.3"=numeric(length(ktests)), "0.4"=numeric(length(ktests)), "0.5"=numeric(length(ktests)), "0.6"=numeric(length(ktests)), "0.7"=numeric(length(ktests)), "0.8"=numeric(length(ktests)))
resultdf_3 = data.frame("K"=ktests,"0.3"=numeric(length(ktests)), "0.4"=numeric(length(ktests)), "0.5"=numeric(length(ktests)), "0.6"=numeric(length(ktests)), "0.7"=numeric(length(ktests)), "0.8"=numeric(length(ktests)))
resultdf_4 = data.frame("K"=ktests,"0.3"=numeric(length(ktests)), "0.4"=numeric(length(ktests)), "0.5"=numeric(length(ktests)), "0.6"=numeric(length(ktests)), "0.7"=numeric(length(ktests)), "0.8"=numeric(length(ktests)))

# setting up empty dataframes to save NMI scores across parameters
CT_NMI = data.frame(matrix(0, ncol = 20, nrow = 68))
SA_NMI = data.frame(matrix(0, ncol = 20, nrow = 68))
FA_NMI = data.frame(matrix(0, ncol = 20, nrow = 46))
vol_NMI = data.frame(matrix(0, ncol = 20, nrow = 14))
clinical_NMI = data.frame(matrix(0, ncol = 20, nrow = 7))

t = 10; 	# Number of Iterations, usually (10~20)

# range of parameters tested
alphas_list <- c("0.3", "0.4", "0.5", "0.6", "0.7", "0.8")

# calculating number of clusters per alpha and K value
# calculating each K for each alpha and saving them

# importing the different data types as individual participant matrices 
subjects <- read.csv("data/Data_for_SNF/ids.csv", header=FALSE)
CT <- read.csv("data/Data_for_SNF/CT_snf.csv", header=FALSE)
volume <- read.csv("data/Data_for_SNF/Volumes_snf.csv", header=FALSE)
FA <- read.csv("data/Data_for_SNF/FA_snf.csv", header=FALSE)
clinical <-  read.csv("data/Data_for_SNF/pond_snf.csv", header=FALSE)

# normalizing measures within each data type using a function from the SNF package
CT = standardNormalization(CT)
volume = standardNormalization(volume)
FA = standardNormalization(FA)
clinical = standardNormalization(clinical)

# creating participant distance matrices using euclidean distances
Dist_CT = dist2(as.matrix(CT),as.matrix(CT));
Dist_volume = dist2(as.matrix(volume),as.matrix(volume));
Dist_FA = dist2(as.matrix(FA),as.matrix(FA));
Dist_clinical = dist2(as.matrix(clinical),as.matrix(clinical));

# for each parameter combination - SNF is used to create a similarity matrix and (1) the optimal number of clusters as well as
# (2) the normalized mutual information (NMI) between each measure and the matrix are calculated and recorded
for(param in 1:length(alphas_list)){
  alpha <- as.numeric(as.character(alphas_list[param]))
  print(alpha)
  for(test in 1:nrow(resultdf_1)){
    K <- resultdf_1$K[test]
    print(test)
    
    # creating participant affinity matrices for each data type
    AM_volume = affinityMatrix(Dist_volume,K, alpha)
    AM_CT = affinityMatrix(Dist_CT,K,alpha)
    AM_FA = affinityMatrix(Dist_FA,K, alpha)
    AM_clinical = affinityMatrix(Dist_clinical,K, alpha)
    
    # Using SNF to create fused similarity matrix
    SNF1 = SNF(list(AM_CT, AM_FA, AM_clinical, AM_volume), K, t)
    
    # calculating optimal number of clusters given the created similarity matrix
    EstClust_1 <-estimateNumberOfClustersGivenGraph(SNF1, NUMC=2:10)[[1]][1] # or [[1]][2]
    EstClust_2 <-estimateNumberOfClustersGivenGraph(SNF1, NUMC=2:10)[[2]][1] # or [[1]][2]
    EstClust_3 <-estimateNumberOfClustersGivenGraph(SNF1, NUMC=2:10)[[3]][1] # or [[1]][2]
    EstClust_4 <-estimateNumberOfClustersGivenGraph(SNF1, NUMC=2:10)[[4]][1] # or [[1]][2]
    
    # adding cluster number to tracking file for the given parameter
    resultdf_1[test,(param +1)]<-EstClust_1
    resultdf_2[test,(param + 1)]<-EstClust_2
    resultdf_3[test,(param +1)]<-EstClust_3
    resultdf_4[test,(param +1)]<-EstClust_4
    
    #calculating NMI scores based on similarity matrix
    SNF1_NMIScores <-rankFeaturesByNMI(list(FA, clinical, volume, CT), SNF1)
    # organizing NMI scores
    FA_NMI[, test] <- SNF1_NMIScores[[1]][1]
    clinical_NMI[, test] <- SNF1_NMIScores[[1]][2]
    vol_NMI[, test] <- SNF1_NMIScores[[1]][3]
    CT_NMI[, test] <- SNF1_NMIScores[[1]][4]
    
    # calculating groups using spectral clustering
    Group <-spectralClustering(SNF1, 4)
    #displayClusters(SNF1,Group)
    
  }
  # setting directory and saving NMI score files
  directory <- (paste("output/1_SNF_analysis/parameter_iterations/alpha_", alpha, sep=""))
  write.csv(CT_NMI, file=file.path(directory, paste("/CT_snf.csv", sep="")))
  write.csv(FA_NMI, file=file.path(directory, paste("/FA_snf.csv", sep="")))
  write.csv(vol_NMI, file=file.path(directory, paste("/vol_snf.csv", sep="")))
  write.csv(clinical_NMI, file=file.path(directory, paste("/clin_snf.csv", sep="")))
  all_scores <- rbind(CT_NMI, FA_NMI)
  all_scores <- rbind(all_scores, vol_NMI)
  all_scores <- rbind(all_scores, clinical_NMI)
  write.csv(all_scores, file=file.path(directory, paste("/all_scores.csv", sep="")))
}

# saving estimated number of clusters files
directory <- ("output/1_SNF_analysis/parameter_iterations/")

write.csv(resultdf_1, file=file.path(directory, paste("Estnumclus_1.csv", sep="")))
write.csv(resultdf_2, file=file.path(directory, paste("Estnumclus_2.csv", sep="")))
write.csv(resultdf_3, file=file.path(directory, paste("Estnumclus_3.csv", sep="")))
write.csv(resultdf_4, file=file.path(directory, paste("Estnumclus_4.csv", sep="")))
