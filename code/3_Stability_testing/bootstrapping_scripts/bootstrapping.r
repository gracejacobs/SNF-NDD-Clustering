####################################################################################
# Script to create bootstrap functions for clustering:
# 1. Create clustering solutions for designated number of permutations
# 2. Calculate the NMI scores for each feature for each permutation
# 3. Determine what the adjusted rand index is between all clustering solutions/the stability of the clusters
# 4. Determine how often each subject is clusted together and the probability that they will be clustered together

# Setting up matrix of which participants will be included in each resampling of 80% of participants
bootstrapping_SNF <- function(numboot, nsub, bootsize){
  

library(SNFtool)
library(cluster)
library(fossil)
library(e1071)

### Bootstrapping results

numboot <- numboot # number of permutations can up to 1000 later
nsub <- nsub #number of subjects
perms <- data.frame(matrix(0, nrow = numboot, ncol = nsub))
bootsize <- bootsize #what percentage of participants do you want to take per permuation
bn <- nsub*bootsize #number of subjects in each bootstrap

# creating output matrices
# agreement <- data.frame(matrix(0, nrow = nsub, ncol = nsub)) # matrix of subject x subject - how often they are clustered together
#num_perms <- data.frame(matrix(0, nrow = nsub, ncol = nsub))
# totagree <- data.frame(matrix(0, nrow = nsub, ncol = nsub)) # how often participants agree with each other
# numinc <- data.frame(matrix(0, nrow = nsub, ncol = nsub)) # number of times each participant is included in a clustering solution
# list_randindex <- data.frame(matrix(NA, nrow = numboot, ncol = numboot)) # creating matrix of adjusted rand indices - coefficient of overlap between each pair of clustering solutions

# 1. Create clustering solutions for designated number of permutations

## making sure none of the permuations have been used before with the same subjects
## creating a matrix of which subjects will be included for each permutation - 1 means they are included, 0 means they are not
print("1. Creating matrix of which participants to include for each permutation")
for(idx in seq(1:numboot)){
  print(idx)
  test <- 0
  while(test == "0"){
    test <- 1
    rnd <- sample(1:nsub, bn) #choosing 80% of subjects
    inc <- data.frame(matrix(0, nrow = 1, ncol = nsub)) # row of all subjects to determine which ones are included in this clustering
    # if the column of inc is not included in rnd, then it will be equal to 0
    for (i in 1:ncol(inc)){
      inc[1,i] <- ifelse(i %in% rnd, 1, 0)
    }
    # checking to see if inc is the same as any of the other perms
    for (row in 1:nrow(perms)){
      bad <- ifelse(perms[row, ] == inc[1, ], 1, 0)
    } 
    if (bad == "1"){ # aka, this is the same clustering solution as any of the previous ones
      test <- 0
    } else {
      test <- 1
    }
  }
  for (i in 1:ncol(inc)){ # adding the row of subjects for the permutation to the permutation matrix
      perms[idx,i] <- inc[1,i]
    }
} 
return(perms)

}

## setting up the insert row function so that I can add rows of 0's for the subjects that were not included in the clustering solution
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

## getting the clustering solutions for all the permuatations using SNF 
clustering <- function(perms, bootsize, K, t, alpha, clusters, CT, FA, clinical, volume){
clus_out <- data.frame(matrix(0, nrow = numboot, ncol = nsub)) # what the clustering is for each permutation
silhouette <- data.frame(matrix(0, nrow = numboot, ncol = nsub))

# Setting SNF parameters
K =K;		# number of neighbors, usually (10~30), usually sample size/10
alpha = alpha;  	# hyperparameter, usually (0.3~0.8)
t = t; 	# Number of Iterations, usually (10~20)
C=clusters #number of clusters you want your participants to group into

# perms <- permutation_matrix
# idx <- 1

print("2. Getting the clustering solutions, as well as silhouette scores for all the permuatations using SNF")
for(idx in 1:numboot){
  print(idx) # permutation number
  subjects <- t(perms[idx, ]) # getting subjects for that permutation
  
  # need to do this for each permutation to re-subset
  temp_FA <- FA
  temp_CT <- CT
  temp_clinical <- clinical
  temp_volume <- volume
  # adding subject number column for later reference
  temp_FA$sub <- subjects
  temp_CT$sub <-subjects
  temp_clinical$sub <-subjects
  temp_volume$sub <-subjects
  # subsetting participants based on permuatations
  temp_FA <- temp_FA[which(temp_FA$sub == "1"), ]
  temp_CT <- temp_CT[which(temp_CT$sub == "1"), ]
  temp_clinical <- temp_clinical[which(temp_clinical$sub == "1"), ]
  temp_volume <- temp_volume[which(temp_volume$sub == "1"), ]
  #removing now unnecessary subject column
  temp_FA$sub <- NULL
  temp_CT$sub <- NULL
  temp_clinical$sub <- NULL
  temp_volume$sub <- NULL
  # SNF clustering with each data type
  temp_CT = standardNormalization(temp_CT)
  temp_volume = standardNormalization(temp_volume)
  temp_FA = standardNormalization(temp_FA)
  temp_clinical = standardNormalization(temp_clinical)
  
  Dist_CT = dist2(as.matrix(temp_CT),as.matrix(temp_CT));
  Dist_volume = dist2(as.matrix(temp_volume),as.matrix(temp_volume));
  Dist_FA = dist2(as.matrix(temp_FA),as.matrix(temp_FA));
  Dist_clinical = dist2(as.matrix(temp_clinical),as.matrix(temp_clinical));

  AM_volume = affinityMatrix(Dist_volume,K, alpha)
  AM_CT = affinityMatrix(Dist_CT,K,alpha) 
  AM_FA = affinityMatrix(Dist_FA,K, alpha)
  AM_clinical = affinityMatrix(Dist_clinical,K, alpha)
  
  SNF1 = SNF(list(AM_CT, AM_volume, AM_FA, AM_clinical), K, t)  
  group = spectralClustering(SNF1, C)
  
  ##getting silhouette plot
  dissim <- 1 - SNF1
  #clusters$groups <- as.integer(group)
  sil <- silhouette(group, dmatrix = dissim)
  
  bn <- as.numeric(nrow(subjects)*0.8)
  silhouette_width <- as.data.frame(matrix(0, ncol = 1, nrow =bn ))
  for (i in 1:nrow(silhouette_width)){
    silhouette_width[i, 1] <- sil[i ,3]
  }
  
  subjects <- as.data.frame(subjects)
  row <- c("0")
  for (i in 1:nrow(subjects)){ #introducing 0 rows to indicate that that subject was not included in clustering
    if (subjects[i,1] == "0"){
      silhouette_width <- insertRow(silhouette_width, row, i)
    } 
  }
  silhouette[idx, ] <- silhouette_width[,1]
  
  ## need to combine with all participants again
  group <- as.data.frame(group)
  subjects <- as.data.frame(subjects)
  row <- c("0")
  for (i in 1:nrow(subjects)){ #introducing 0 rows to indicate that that subject was not included in clustering
    if (subjects[i,1] == "0"){
      group <- insertRow(group, row, i)
    } 
  }
  clus_out[idx, ] <- group[,1] # adding the clustering solution to the final output matrix clus_out
  
}
  clus_sil <- rbind(clus_out, silhouette)
  return(clus_sil)
}

# 2. Calculate the NMI scores for each feature for each permutation

# getting NMI scores for each of the cluster permutations
NMI_scores <- function(perms, bootsize, K, t, alpha, clusters, CT, FA, clinical, volume){
  clus_out <- data.frame(matrix(0, nrow = numboot, ncol = nsub)) # what the clustering is for each permutation
  
  # Setting SNF parameters
  K =K;		# number of neighbors, usually (10~30), usually sample size/10
  alpha = alpha;  	# hyperparameter, usually (0.3~0.8)
  t = t; 	# Number of Iterations, usually (10~20)
  C=clusters #number of clusters you want your participants to group into

  CT_NMI <- data.frame(matrix(0, nrow = length(CT), ncol = numboot))
  clin_NMI <- data.frame(matrix(0, nrow = length(clinical), ncol = numboot))
  FA_NMI <- data.frame(matrix(0, nrow = length(FA), ncol = numboot))
  vol_NMI <- data.frame(matrix(0, nrow = length(volume), ncol = numboot))
  
print("2. Getting the clustering solutions for all the permuatations using SNF, and then getting the NMI scores")
for(idx in 1:numboot){
  print(idx) # permutation number
  subjects <- t(perms[idx, ]) # getting subjects for that permutation
  
  # need to do this each permutation to re-subset
  temp_FA <- FA
  temp_CT <- CT
  temp_clinical <- clinical
  temp_volume <- volume
  # adding subject number column for later reference
  temp_FA$sub <- subjects
  temp_CT$sub <-subjects
  temp_clinical$sub <-subjects
  temp_volume$sub <-subjects
  # subsetting participants based on permuatations
  temp_FA <- temp_FA[which(temp_FA$sub == "1"), ]
  temp_CT <- temp_CT[which(temp_CT$sub == "1"), ]
  temp_clinical <- temp_clinical[which(temp_clinical$sub == "1"), ]
  temp_volume <- temp_volume[which(temp_volume$sub == "1"), ]
  #removing now unnecessary subject column
  temp_FA$sub <- NULL
  temp_CT$sub <- NULL
  temp_clinical$sub <- NULL
  temp_volume$sub <- NULL
  # SNF clustering with each data type
  temp_CT = standardNormalization(temp_CT)
  temp_volume = standardNormalization(temp_volume)
  temp_FA = standardNormalization(temp_FA)
  temp_clinical = standardNormalization(temp_clinical)
  
  Dist_CT = dist2(as.matrix(temp_CT),as.matrix(temp_CT));
  Dist_volume = dist2(as.matrix(temp_volume),as.matrix(temp_volume));
  Dist_FA = dist2(as.matrix(temp_FA),as.matrix(temp_FA));
  Dist_clinical = dist2(as.matrix(temp_clinical),as.matrix(temp_clinical));
  
  AM_volume = affinityMatrix(Dist_volume,K, alpha)
  AM_CT = affinityMatrix(Dist_CT,K,alpha) 
  AM_FA = affinityMatrix(Dist_FA,K, alpha)
  AM_clinical = affinityMatrix(Dist_clinical,K, alpha)
  
  SNF1 = SNF(list(AM_CT, AM_volume, AM_FA, AM_clinical), K, t)  
  
  ## getting NMI scores 
  SNF1_NMIScores <-rankFeaturesByNMI(list(temp_FA, temp_clinical, temp_volume, temp_CT), SNF1)
  
  CT_NMI[,idx] <- SNF1_NMIScores[[1]][4]
  clin_NMI[,idx] <- SNF1_NMIScores[[1]][2]
  FA_NMI[,idx] <- SNF1_NMIScores[[1]][1]
  vol_NMI[,idx] <- SNF1_NMIScores[[1]][3]
  
}
## need to combine NMIs first
All_NMI_scores <- rbind(clin_NMI, vol_NMI)
All_NMI_scores <- rbind(All_NMI_scores, FA_NMI)
All_NMI_scores <- rbind(All_NMI_scores, CT_NMI)

return(All_NMI_scores)
}

# 3. Determine what the adjusted rand index is between all clustering solutions/the stability of the clusters

## CHecking the adjusted rand index (coefficient of how similar the clustering is) for each pair of clustering solutions
# 0 means no coherance between clustering up to 1

# for each clustering solution, compare to all other clustering solutions

stability <- function(clus_out, perms){
print("3. Comparing each clustering solution to get the adjusted rand index")
list_randindex <- data.frame(matrix(NA, nrow = numboot, ncol = numboot)) # creating matrix of adjusted rand indices - coefficient of overlap between each pair of clustering solutions
list_adjrandindex <- data.frame(matrix(NA, nrow = numboot, ncol = numboot)) # creating matrix of adjusted rand indices - coefficient of overlap between each pair of clustering solutions
  
  
for (rep in 1:(numboot - 1)){
  print(rep)
  for (idx in (rep + 1):numboot) {
    subs_list <- perms[c(rep, idx), ] # getting list of subjects for that permutation
    subs_list[3, ] <- seq(1:ncol(subs_list)) # creating a column of subject number
    # finding the number of participants that are a part of both solutions
    subs_list <- subs_list[, which(subs_list[1, ] != "0" & subs_list[2, ] != "0")] 
    subs_list <-as.numeric(t(subs_list[3, ]))
    # getting clusters for those overlapping subjects
    c1 <- clus_out[rep, c(subs_list)]
    c2 <- clus_out[idx, c(subs_list)]
    c1 <- as.numeric(c1)
    c2 <- as.numeric(c2)
    # getting the adjusted rand index of the two clustering lists
    list_randindex[rep, idx] <- rand.index(t(c1), t(c2))
    list_adjrandindex[rep, idx] <- adj.rand.index(t(c1), t(c2))
  }
}

list_randindex <- cbind(list_randindex, list_adjrandindex)
return(list_randindex)

}

# 4. Determine how often each subject is clusted together and the probability that they will be clustered together

## get the agreement matrix - how often two subjects are clustered together
## for each permutation, for each subject

percent_agree <- function(clus_out){
print("4. Determine how often each subject is clustered together")
  
  agreement <- data.frame(matrix(0, nrow = nsub, ncol = nsub))
  totagree <- data.frame(matrix(0, nrow = nsub, ncol = nsub)) # how often participants agree with each other
  numinc <- data.frame(matrix(0, nrow = nsub, ncol = nsub)) # number of times each participant is included in a clustering solution
  
for (perm in 1:numboot){
  print(perm) # permutation number
  data <- as.data.frame(clus_out[perm, ])
  for (idx in 1:nsub){
    if (data[ ,idx] > 0){ # if they are included in the clustering
      matched <- as.data.frame(t(data[1, ])) # permutation list of subjects
      matched$num <- seq(1:nrow(matched)) # creating a column of subject number
      comp_matched <- matched[which(matched[ ,1] == matched[idx, 1]), ] # getting which subjects have the same cluster number of that subject
      numlist <- as.numeric(t(comp_matched$num)) # creating list of them
      totagree[idx, c(numlist)] <- totagree[idx, c(numlist)] + 1 # increasing number in the totagree matrix if they do match
      
      inc <- matched[which(matched[ ,1] != "0"), ] # finding out which subjects were included in clustering for that perm with that subject
      inclist <- as.numeric(t(inc$num))
      numinc[idx, c(inclist)] <- numinc[idx, c(inclist)] + 1 ## increasing number in the numinc matrix if they do match
    }
  }
}

# percent of time each subject is clustered together
percent_agree <- totagree/numinc
print("Done bootstrapping 4/4 steps")

return(percent_agree) 
}
