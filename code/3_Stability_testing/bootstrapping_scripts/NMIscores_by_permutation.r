# getting NMI scores for each permutation

NMI_scores <- function(perms, bootsize, K, t, alpha, clusters, CT, FA, clinical, volume){
  clus_out <- data.frame(matrix(0, nrow = numboot, ncol = nsub)) # what the clustering is for each permutation
  
  # Setting SNF parameters
  K =K;		# number of neighbors, usually (10~30), usually sample size/10
  alpha = alpha;  	# hyperparameter, usually (0.3~0.8)
  t = t; 	# Number of Iterations, usually (10~20)
  C=clusters #number of clusters you want your participants to group into
  
  insertRow <- function(existingDF, newrow, r) {
    existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
    existingDF[r,] <- newrow
    existingDF
  }
  
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
