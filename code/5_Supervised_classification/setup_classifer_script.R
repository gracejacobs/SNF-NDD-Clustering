### Creates function for random forest classifier
## divides participants into an 80/20 train test split
library(e1071)
library(rpart)
library(fossil)
library(caret)

# function to divide the data into training and testing samples across permutations
dividing_data <- function(numboot, nsub, bootsize){
  
  perms <- data.frame(matrix(0, nrow = numboot, ncol = nsub))
  bn <- nsub*bootsize #number of subjects in each bootstrap

  train <- data.frame(matrix(0, nrow = numboot, ncol = nsub))
  testing <- data.frame(matrix(0, nrow = numboot, ncol = nsub))
  testing_labels <- data.frame(matrix(0, nrow = numboot, ncol = nsub))

  print("1. Creating matrices of which participants to include for training and testing (80%)")
  for(idx in seq(1:numboot)){
    print(idx)
    test <- 0
    while(test == "0"){
      test <- 1
      rnd_1 <- sample(1:nsub, bn)
      
      inc <- data.frame(matrix(0, nrow = 1, ncol = nsub)) # row of all subjects to determine which ones are included in this clustering
      # if the column of inc is not included in rnd, then it will be equal to 0
      for (i in 1:ncol(inc)){
        #inc[1,i] <- ifelse(i %in% rnd_1 | i %in% rnd_2 | i %in% rnd_3 | i %in% rnd_4, 1, 0)
        inc[1,i] <- ifelse(i %in% rnd_1, 1, 0)
      }
      
      # checking to see if inc is the same as any of the other perms
      bad <- 0 
      for (row in 1:nrow(train)){
        bad <- ifelse(identical(perms[row, ], inc[1, ]), bad + 1, bad)
      } 
      
      if (bad > 0){ # aka, this is the same clustering solution as any of the previous ones
        test <- 0
      } else {
        test <- 1
      }
    }
    
    for (i in 1:ncol(inc)){ # adding the row of subjects for the permutation to the permutation matrix
      train[idx,i] <- inc[1,i]
    }
    
    ## Setting up testing matrix
    testing_peeps <- data.frame(matrix(0, nrow = 1, ncol = nsub))
    for (i in 1:ncol(inc)){
      testing_peeps[1,i] <- ifelse(inc[1,i] == "0", 1, 0)
    }
    for (i in 1:ncol(testing_peeps)){ # adding the row of subjects for the permutation to the permutation matrix
      testing[idx,i] <- testing_peeps[1,i]
    }
  } 
  # putting both together so that I can return them
  both <- rbind(train, testing)
  return(both)
}
################################################################################################################################
################################################################################################################################
