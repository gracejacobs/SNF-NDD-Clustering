### Running SNF on the POND data
## This script is run after Running_parameter_iterations.r in order to determine the optimal 
# cluster number, hyperparameter (alpha), and nearest neighbors parameter (K)

library(SNFtool)
library(psych)
library(dunn.test)
library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)
library(fossil)
library(corrplot)
library(cluster)
library(MASS)

# source code for function to determine data integration and clustering across resampling using SNF and spectral clustering
source("code/1_Running_SNF_analysis/robust_core-clustering_function.R")

# set up for SNF analysis
# importing the different data types as individual participant matrices 
subjects <- read.csv("data/Data_for_SNF/ids.csv", header=FALSE) #list of participants
CT <- read.csv("data/Data_for_SNF/CT_snf.csv", header=FALSE) # cortical thickness (n=68)
volume <- read.csv("data/Data_for_SNF/Volumes_snf.csv", header=FALSE) # subcortical volumes (n=14)
FA <- read.csv("data/Data_for_SNF/FA_snf.csv", header=FALSE) # white matter fractional anisotropy (n=46)
clinical <-  read.csv("data/Data_for_SNF/pond_snf.csv", header=FALSE) # behavioural measures (n=7)

# normalizing measures within each data type using a function from the SNF package
CT = standardNormalization(CT)
volume = standardNormalization(volume)
FA = standardNormalization(FA)
clinical = standardNormalization(clinical)

# setting the parameters (finalized after comparisons using Running_parameter_iterations.r )
K =18;		# number of neighbors, usually (10~30), usually sample size/10
alpha = 0.8;  	# hyperparameter, usually (0.3~0.8)
t = 10; 	# Number of Iterations, usually (10~20) 

# creating participant distance matrices using euclidean distances
Dist_CT = dist2(as.matrix(CT),as.matrix(CT));
Dist_volume = dist2(as.matrix(volume),as.matrix(volume));
Dist_FA = dist2(as.matrix(FA),as.matrix(FA));
Dist_clinical = dist2(as.matrix(clinical),as.matrix(clinical));

# creating participant affinity matrices within each data type
AM_volume = affinityMatrix(Dist_volume,K, alpha)
AM_CT = affinityMatrix(Dist_CT,K,alpha) 
AM_FA = affinityMatrix(Dist_FA,K, alpha)
AM_clinical = affinityMatrix(Dist_clinical,K, alpha)

### Calculating similarity matrix and spectral clustering groups
# setting output directory
directory <- ("output/1_SNF_analysis/")

C = 4 #setting cluster number
# calling function to integrate data types using SNF and cluster participants using spectral clustering 
#across resampling 80% of participants 1000 times
robust.W = RobustCoreClusteringMatrix(feature.affinity.mat.list = list(AM_FA, AM_CT, AM_volume, AM_clinical),
                                      exp.num.samples = 1000, num.clusts = C)
#Two matrices - Dense Core Cluster Matrix and Sparse Core Cluster Matrix
dense <- robust.W[1]
dense <- matrix(unlist(dense), ncol = 176, byrow = TRUE)
sparse <- robust.W[2]
sparse <- matrix(unlist(sparse), ncol = 176, byrow = TRUE)

# displaying clusters
displayClustersWithHeatmap(dense, spectralClustering(dense, C))
displayClustersWithHeatmap(sparse, spectralClustering(sparse, C))

# saving an image of the cluster heatmap
png('SN_dense_0.8_18_1000perms.png')
displayClustersWithHeatmap(dense, spectralClustering(dense, C))
dev.off()

# calculating normalized mutual information (NMI) based off of the original data types and the clustering similarity matrix
SNF_NMIScores <-rankFeaturesByNMI(list(volume, FA, CT, clinical), dense) 
SNF_NMIScores 

# separating and organizing scores
vol_scores <- as.data.frame(SNF_NMIScores[[1]][1])
FA_scores <- as.data.frame(SNF_NMIScores[[1]][2])
CT_scores <- as.data.frame(SNF_NMIScores[[1]][3])
clin_scores <- as.data.frame(SNF_NMIScores[[1]][4])
names(FA_scores) <- c("NMI")
names(clin_scores) <- c("NMI")
names(vol_scores) <- c("NMI")
names(CT_scores) <- c("NMI")

all_scores <- rbind(CT_scores, FA_scores)
all_scores <- rbind(all_scores, vol_scores)
all_scores <- rbind(all_scores, clin_scores)

#saving csv of NMI scores and clustering similarity matrix
write.csv(all_scores, file=file.path(directory, paste("all_scores_k18_0.8_1000perms.csv", sep="")))
write.matrix(dense, file=file.path(directory, paste(name, sep="")))

# Find cluster labels of individuals using the robust clustering similarity matrix
robust.groups.df = RobustCoreClusteringClusters(core.clustering.list = robust.W,num.clusts = C,verbose = T)
clusters <- cbind(subjects, robust.groups.df)
table(clusters$groups)

## calculating silouette width for each participant and the silhouette plot
dissim <- 1 - dense
dissim <- as.matrix(dissim)

# reordering clusters from generally least impaired to most impaired (this was determined later when comparing clusters)
#clusters <- read.csv("output/1_SNF_analysis/4clust_groups_k18_0.8_1000perms.csv")
clusters$groups_2 <- ifelse(clusters$groups == "1", "2", NA)
clusters$groups_2 <- ifelse(clusters$groups == "2", "4", clusters$groups_2)
clusters$groups_2 <- ifelse(clusters$groups == "3", "3", clusters$groups_2)
clusters$groups_2 <- ifelse(clusters$groups == "4", "1", clusters$groups_2)
clusters$groups <- clusters$groups_2
clusters$groups <- as.integer(clusters$groups)

sil <- silhouette(clusters$groups, dmatrix = dissim)

# saving the silhouette plot
name=("Silhouette_plot_0.8_18.png")
png(name)
plot(sil, col = c("palevioletred", "darkkhaki", "turquoise3", "slateblue"))
dev.off()

# adding silhouette widths for each individual
clusters$silhouette_width <- 0
for (i in 1:176){
  clusters[i, 4] <- sil[i ,3]
}

write.csv(clusters, file=file.path(directory, paste("4clust_groups_k18_0.8_1000perms.csv", sep="")))


# Check if anyone doesn't reliably cluster 
unreliable.patients = UnrelClustPatientsFullCohort(core.clustering.list = robust.W,verbose = TRUE)

# Check if anyone doesn't reliably cluster within their given cluster group
unreliable.patients.by.grp = UnrelClustPatientsByGrp(core.clustering.list = robust.W,id.group.df = robust.groups.df,verbose = TRUE)




