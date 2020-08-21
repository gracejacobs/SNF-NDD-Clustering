
## Getting the stability of each clustering solution

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


#### getting hamming distance for each permutation
hamming <- function(clus_out, perms){
  print("3. Comparing each clustering solution to get the hamming distance")
  hamming_distance <- data.frame(matrix(NA, nrow = numboot, ncol = numboot)) # creating matrix of adjusted rand indices - coefficient of overlap between each pair of clustering solutions
  
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
      hamming_distance[rep, idx] <- hamming.distance(c1, c2)
    }
  }
  return(hamming_distance)
}