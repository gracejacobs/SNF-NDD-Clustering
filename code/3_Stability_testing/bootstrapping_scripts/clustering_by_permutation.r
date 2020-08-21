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
  
  ## setting up the insert row function so that I can add rows of 0's for the subjects that were not included in the clustering solution
  insertRow <- function(existingDF, newrow, r) {
    existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
    existingDF[r,] <- newrow
    existingDF
  }
  
  print("2. Getting the clustering solutions, as well as silhouette scores for all the permuatations using SNF")
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
    group = spectralClustering(SNF1, C)
    
    ##getting silhouette
    
    dissim <- 1 - SNF1
    #clusters$groups <- as.integer(group)
    sil <- silhouette(group, dmatrix = dissim)
    # setwd("/projects/gjacobs/POND/Running_SNF/SNF_Dec12_fullqc/Bootstrapping_k23_0.8/Silhouette_plots")
    # png('4_clusters_dense_0.8_23.png')
    # plot(sil, col = c("red", "green", "blue", "purple"))
    # dev.off()
    
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
