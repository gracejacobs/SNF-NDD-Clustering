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
  percent_agree <<- totagree/numinc
  print("Done bootstrapping 4/4 steps")
  
  return(percent_agree) 
}