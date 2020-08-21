## creating permuataion matrix for bootstrapping analyses

bootstrapping_SNF <- function(numboot, nsub, bootsize){
  
  numboot <- numboot # number of permutations can up to 1000 later
  nsub <- nsub #number of subjects
  perms <- data.frame(matrix(0, nrow = numboot, ncol = nsub))
  bootsize <- bootsize #what percentage of participants do you want to take per permuation
  bn <- nsub*bootsize #number of subjects in each bootstrap
  
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