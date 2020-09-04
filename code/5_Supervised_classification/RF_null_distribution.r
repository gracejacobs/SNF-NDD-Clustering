## Creating a null distribution of classification model metrics to determine performance significance

library(randomForest)
library(caret)
library(data.table)

## importing cluster labels
labels <- read.csv("output/1_SNF_analysis/4clust_groups_k18_0.8_1000perms.csv")
labels$groups_2 <- ifelse(labels$groups == "1", "2", NA)
labels$groups_2 <- ifelse(labels$groups == "2", "4", labels$groups_2)
labels$groups_2 <- ifelse(labels$groups == "3", "3", labels$groups_2)
labels$groups_2 <- ifelse(labels$groups == "4", "1", labels$groups_2)
labels$groups <- labels$groups_2
labels <- as.data.frame(labels$groups)

ids <- read.csv("data/Data_with_headers/ids.csv") # participant subject identifiers
measures <- read.csv("data/Data_with_headers/all_measures.csv")
measures <- measures[, c(2:129, 198:204)]

# setting up the classifier
source("code/5_Supervised_classification/setup_classifer_script.R")

numboot <- 500 # number of permutations can up to 1000 later
nsub <- 176 #number of subjects
bootsize <- 0.8 #what percentage of participants do you want to take per permuation

## getting the testing and training set
both <- dividing_data(numboot=numboot, nsub=nsub, bootsize=bootsize, labels=labels, measures=measures, ids=ids)
train <- both[1:numboot, ]
test <- both[(1+numboot):(2*numboot), ]

###########################################################################################
# All features included in the classification model

names <- read.csv("output/1_SNF_analysis/all_scores_k18_0.8_1000perms.csv", header = TRUE)
names <- names[1:135, 3]
#adjusting tract names
names <- gsub('-', '.', names)
names <- gsub('/', '.', names)
names <- as.data.frame(names)

measures <- measures[,colnames(measures) %in% names] 

# adding the cluster labels to the measures
lm <- cbind(labels, measures)
# adding the cluster labels to the ids
lm <- cbind(ids, lm)
# ordering measures by cluster labels
lm <- lm[order(lm$`labels$groups`), ]
lm$x <- NULL

perms <- data.frame(matrix(0, nrow = numboot, ncol = nsub))
bn <- nsub*bootsize #number of subjects in each bootstrap

##creating matrices to keep track of accuracies
acc_overall <- data.frame(matrix(0, nrow = 7, ncol = numboot))
acc_class <- data.frame(matrix(0, nrow = 1, ncol = 11))
list_adjrandindex <- data.frame(matrix(NA, nrow = numboot, ncol = 1)) # creating matrix of adjusted rand indices - coefficient of overlap between each pair of clustering solutions

### Running the actual random forest classifer
for(idx in 1:numboot){
  print(idx) # permutation number
  subjects <- t(train[idx, ]) # getting subjects for that permutation
  
  training_set <- cbind(subjects, lm)
  training_set <- training_set[which(training_set[,1] == "1"),] ## only getting the people for this perm
  ## shuffling labels
  training_set[,2] <- sample(training_set$`labels$groups`, replace = FALSE)
  labels_train <- as.factor(training_set[ ,2])
  training_set[,1:2] <- NULL
  
  ## setting up testing dataset - don't want to shuffle these
  subjects <- t(test[idx, ]) # getting subjects for that permutation
  testing_set <- cbind(subjects, lm)
  testing_set <- testing_set[which(testing_set[,1] == "1"),]
  labels_test <- as.factor(testing_set[ ,2])
  testing_set[,1:2] <- NULL
  # running the model
  model.rf = randomForest(labels_train ~. , data=training_set, ntree=501, mtry=80, importance=TRUE)

  rfpredict <- predict(model.rf, testing_set)
  stats <- confusionMatrix(rfpredict, labels_test) 
  
  acc_bygroup_perm <- as.data.frame(stats[4][1])
  names(acc_class) <- names(acc_bygroup_perm)
  acc_class <- rbind(acc_class, acc_bygroup_perm)
  
  row.names(acc_overall) <- row.names(as.data.frame(stats[3][1]))
  acc_overall[ , idx] <- as.data.frame(stats[3][1])
  list_adjrandindex[idx, 1] <- adj.rand.index(rfpredict, labels_test)
}


# calculating a significance cut-off (p<0.05) for adjusted rand index 
names(list_adjrandindex) <- c("index")
adjrandindex <- sort(list_adjrandindex$index)
quantile(adjrandindex, 0.95)
# calculating a significant (p<0.05) overall model accuracy
acc_overall <- as.data.frame(t(acc_overall))
acc_all <- acc_overall
acc_all <- sort(acc_all$Accuracy)
quantile(acc_all, 0.95)
# calculating a significance cut-off (p<0.05) for sensitivity and specificity for each group
acc_class <-acc_class[-1,] 
acc_class$label <- rep(1:4, times=100, each=1)
acc_class$label <- as.factor(acc_class$label)
acc_1 <- acc_class[which(acc_class$label == "1"), ]
acc_2 <- acc_class[which(acc_class$label == "2"), ]
acc_3 <- acc_class[which(acc_class$label == "3"), ]
acc_4 <- acc_class[which(acc_class$label == "4"), ]
# sensitivity
acc_1 <- acc_1[order(acc_1$byClass.Sensitivity), ]
acc_2 <- acc_2[order(acc_2$byClass.Sensitivity), ]
acc_3 <- acc_3[order(acc_3$byClass.Sensitivity), ]
acc_4 <- acc_4[order(acc_4$byClass.Sensitivity), ]
quantile(acc_1$byClass.Sensitivity, 0.95)
quantile(acc_2$byClass.Sensitivity, 0.95)
quantile(acc_3$byClass.Sensitivity, 0.95)
quantile(acc_4$byClass.Sensitivity, 0.95)
# specificity
acc_1 <- acc_class[which(acc_class$label == "1"), ]
acc_2 <- acc_class[which(acc_class$label == "2"), ]
acc_3 <- acc_class[which(acc_class$label == "3"), ]
acc_4 <- acc_class[which(acc_class$label == "4"), ]
acc_1 <- acc_1[order(acc_1$byClass.Specificity), ]
acc_2 <- acc_2[order(acc_2$byClass.Specificity), ]
acc_3 <- acc_3[order(acc_3$byClass.Specificity), ]
acc_4 <- acc_4[order(acc_4$byClass.Specificity), ]
quantile(acc_1$byClass.Specificity, 0.95)
quantile(acc_2$byClass.Specificity, 0.95)
quantile(acc_3$byClass.Specificity, 0.95)
quantile(acc_4$byClass.Specificity, 0.95)

###########################################################################################
# including top 2 measures from each
names <- read.csv("output/1_SNF_analysis/all_scores_k18_0.8_1000perms.csv", header = TRUE)
names <- names[1:135, 3]
#adjusting tract names
names <- gsub('-', '.', names)
names <- gsub('/', '.', names)
#setdiff(colnames(measures), names)
names <- as.data.frame(names)

names <- names[which(names$names == "ADHD_I_SUB" | names$names == "ADHD_HI_SUB"| names$names == "R_insula_thickavg"| names$names == "R_parstriangularis_thickavg"
                     | names$names =="Rpal"| names$names =="Rput"| names$names == "ALIC.L"| names$names == "RLIC.R"),]

measures <- measures[,colnames(measures) %in% names] 
names(measures)
# checking measure variance
nzv <- nearZeroVar(measures, saveMetrics = TRUE)
sum(nzv$nzv == TRUE)

# adding the cluster labels to the measures
lm <- cbind(labels, measures)
# adding the cluster labels to the ids
lm <- cbind(ids, lm)
# ordering measures by cluster labels
lm <- lm[order(lm$`labels$groups`), ]
lm$x <- NULL

perms <- data.frame(matrix(0, nrow = numboot, ncol = nsub))
bn <- nsub*bootsize #number of subjects in each bootstrap

##creating matrices to keep track of accuracies
acc_overall <- data.frame(matrix(0, nrow = 7, ncol = numboot))
acc_class <- data.frame(matrix(0, nrow = 1, ncol = 11))
list_adjrandindex <- data.frame(matrix(NA, nrow = numboot, ncol = 1)) 

### Running the actual random forest classifer
for(idx in 1:numboot){
  print(idx) # permutation number
  subjects <- t(train[idx, ]) # getting subjects for that permutation
  
  training_set <- cbind(subjects, lm)
  training_set <- training_set[which(training_set[,1] == "1"),] ## only getting the people for this perm
  ## shuffling labels
  training_set[,2] <- sample(training_set$`labels$groups`, replace = FALSE)
  labels_train <- as.factor(training_set[ ,2])
  training_set[,1:2] <- NULL
  
  ## setting up testing dataset - don't want to shuffle these
  subjects <- t(test[idx, ]) # getting subjects for that permutation
  testing_set <- cbind(subjects, lm)
  testing_set <- testing_set[which(testing_set[,1] == "1"),]
  labels_test <- as.factor(testing_set[ ,2])
  testing_set[,1:2] <- NULL
  # running the model
  model.rf = randomForest(labels_train ~. , data=training_set, ntree=501, mtry=4, importance=TRUE)
  
  rfpredict <- predict(model.rf, testing_set)
  stats <- confusionMatrix(rfpredict, labels_test) 
  
  acc_bygroup_perm <- as.data.frame(stats[4][1])
  names(acc_class) <- names(acc_bygroup_perm)
  acc_class <- rbind(acc_class, acc_bygroup_perm)
  
  row.names(acc_overall) <- row.names(as.data.frame(stats[3][1]))
  acc_overall[ , idx] <- as.data.frame(stats[3][1])
  list_adjrandindex[idx, 1] <- adj.rand.index(rfpredict, labels_test)
}


# calculating a significance cut-off (p<0.05) for adjusted rand index 
names(list_adjrandindex) <- c("index")
adjrandindex <- sort(list_adjrandindex$index)
quantile(adjrandindex, 0.95)
# calculating a significant (p<0.05) overall model accuracy
acc_overall <- as.data.frame(t(acc_overall))
acc_all <- acc_overall
acc_all <- sort(acc_all$Accuracy)
quantile(acc_all, 0.95)
# calculating a significance cut-off (p<0.05) for sensitivity and specificity for each group
acc_class <-acc_class[-1,] 
acc_class$label <- rep(1:4, times=100, each=1)
acc_class$label <- as.factor(acc_class$label)
acc_1 <- acc_class[which(acc_class$label == "1"), ]
acc_2 <- acc_class[which(acc_class$label == "2"), ]
acc_3 <- acc_class[which(acc_class$label == "3"), ]
acc_4 <- acc_class[which(acc_class$label == "4"), ]
# sensitivity
acc_1 <- acc_1[order(acc_1$byClass.Sensitivity), ]
acc_2 <- acc_2[order(acc_2$byClass.Sensitivity), ]
acc_3 <- acc_3[order(acc_3$byClass.Sensitivity), ]
acc_4 <- acc_4[order(acc_4$byClass.Sensitivity), ]
quantile(acc_1$byClass.Sensitivity, 0.95)
quantile(acc_2$byClass.Sensitivity, 0.95)
quantile(acc_3$byClass.Sensitivity, 0.95)
quantile(acc_4$byClass.Sensitivity, 0.95)
# specificity
acc_1 <- acc_class[which(acc_class$label == "1"), ]
acc_2 <- acc_class[which(acc_class$label == "2"), ]
acc_3 <- acc_class[which(acc_class$label == "3"), ]
acc_4 <- acc_class[which(acc_class$label == "4"), ]
acc_1 <- acc_1[order(acc_1$byClass.Specificity), ]
acc_2 <- acc_2[order(acc_2$byClass.Specificity), ]
acc_3 <- acc_3[order(acc_3$byClass.Specificity), ]
acc_4 <- acc_4[order(acc_4$byClass.Specificity), ]
quantile(acc_1$byClass.Specificity, 0.95)
quantile(acc_2$byClass.Specificity, 0.95)
quantile(acc_3$byClass.Specificity, 0.95)
quantile(acc_4$byClass.Specificity, 0.95)

###########################################################################################
# including top 35 features
names <- read.csv("output/1_SNF_analysis/all_scores_k18_0.8_1000perms.csv", header = TRUE)
names <- names[1:135, 3]
#adjusting tract names
names <- gsub('-', '.', names)
names <- gsub('/', '.', names)
#setdiff(colnames(measures), names)
names <- as.data.frame(names)

names <- names[1:35, ]

measures <- measures[,colnames(measures) %in% names] 
names(measures)
# checking measure variance
nzv <- nearZeroVar(measures, saveMetrics = TRUE)
sum(nzv$nzv == TRUE)

# adding the cluster labels to the measures
lm <- cbind(labels, measures)
# adding the cluster labels to the ids
lm <- cbind(ids, lm)
# ordering measures by cluster labels
lm <- lm[order(lm$`labels$groups`), ]
lm$x <- NULL

perms <- data.frame(matrix(0, nrow = numboot, ncol = nsub))
bn <- nsub*bootsize #number of subjects in each bootstrap

##creating matrices to keep track of accuracies
acc_overall <- data.frame(matrix(0, nrow = 7, ncol = numboot))
acc_class <- data.frame(matrix(0, nrow = 1, ncol = 11))
list_adjrandindex <- data.frame(matrix(NA, nrow = numboot, ncol = 1)) 

### Running the actual random forest classifer
for(idx in 1:numboot){
  print(idx) # permutation number
  subjects <- t(train[idx, ]) # getting subjects for that permutation
  
  training_set <- cbind(subjects, lm)
  training_set <- training_set[which(training_set[,1] == "1"),] ## only getting the people for this perm
  ## shuffling labels
  training_set[,2] <- sample(training_set$`labels$groups`, replace = FALSE)
  labels_train <- as.factor(training_set[ ,2])
  training_set[,1:2] <- NULL
  
  ## setting up testing dataset - don't want to shuffle these
  subjects <- t(test[idx, ]) # getting subjects for that permutation
  testing_set <- cbind(subjects, lm)
  testing_set <- testing_set[which(testing_set[,1] == "1"),]
  labels_test <- as.factor(testing_set[ ,2])
  testing_set[,1:2] <- NULL
  # running the model
  model.rf = randomForest(labels_train ~. , data=training_set, ntree=501, mtry=15, importance=TRUE)
  
  rfpredict <- predict(model.rf, testing_set)
  stats <- confusionMatrix(rfpredict, labels_test) 
  
  acc_bygroup_perm <- as.data.frame(stats[4][1])
  names(acc_class) <- names(acc_bygroup_perm)
  acc_class <- rbind(acc_class, acc_bygroup_perm)
  
  row.names(acc_overall) <- row.names(as.data.frame(stats[3][1]))
  acc_overall[ , idx] <- as.data.frame(stats[3][1])
  list_adjrandindex[idx, 1] <- adj.rand.index(rfpredict, labels_test)
}

# calculating a significance cut-off (p<0.05) for adjusted rand index 
names(list_adjrandindex) <- c("index")
adjrandindex <- sort(list_adjrandindex$index)
quantile(adjrandindex, 0.95)
# calculating a significant (p<0.05) overall model accuracy
acc_overall <- as.data.frame(t(acc_overall))
acc_all <- acc_overall
acc_all <- sort(acc_all$Accuracy)
quantile(acc_all, 0.95)
# calculating a significance cut-off (p<0.05) for sensitivity and specificity for each group
acc_class <-acc_class[-1,] 
acc_class$label <- rep(1:4, times=100, each=1)
acc_class$label <- as.factor(acc_class$label)
acc_1 <- acc_class[which(acc_class$label == "1"), ]
acc_2 <- acc_class[which(acc_class$label == "2"), ]
acc_3 <- acc_class[which(acc_class$label == "3"), ]
acc_4 <- acc_class[which(acc_class$label == "4"), ]
# sensitivity
acc_1 <- acc_1[order(acc_1$byClass.Sensitivity), ]
acc_2 <- acc_2[order(acc_2$byClass.Sensitivity), ]
acc_3 <- acc_3[order(acc_3$byClass.Sensitivity), ]
acc_4 <- acc_4[order(acc_4$byClass.Sensitivity), ]
quantile(acc_1$byClass.Sensitivity, 0.95)
quantile(acc_2$byClass.Sensitivity, 0.95)
quantile(acc_3$byClass.Sensitivity, 0.95)
quantile(acc_4$byClass.Sensitivity, 0.95)
# specificity
acc_1 <- acc_class[which(acc_class$label == "1"), ]
acc_2 <- acc_class[which(acc_class$label == "2"), ]
acc_3 <- acc_class[which(acc_class$label == "3"), ]
acc_4 <- acc_class[which(acc_class$label == "4"), ]
acc_1 <- acc_1[order(acc_1$byClass.Specificity), ]
acc_2 <- acc_2[order(acc_2$byClass.Specificity), ]
acc_3 <- acc_3[order(acc_3$byClass.Specificity), ]
acc_4 <- acc_4[order(acc_4$byClass.Specificity), ]
quantile(acc_1$byClass.Specificity, 0.95)
quantile(acc_2$byClass.Specificity, 0.95)
quantile(acc_3$byClass.Specificity, 0.95)
quantile(acc_4$byClass.Specificity, 0.95)

###########################################################################################

# including top 10 features
names <- read.csv("output/1_SNF_analysis/all_scores_k18_0.8_1000perms.csv", header = TRUE)
names <- names[1:135, 3]
#adjusting tract names
names <- gsub('-', '.', names)
names <- gsub('/', '.', names)
#setdiff(colnames(measures), names)
names <- as.data.frame(names)

names <- names[1:10, ]

measures <- measures[,colnames(measures) %in% names] 
names(measures)
# checking measure variance
nzv <- nearZeroVar(measures, saveMetrics = TRUE)
sum(nzv$nzv == TRUE)

# adding the cluster labels to the measures
lm <- cbind(labels, measures)
# adding the cluster labels to the ids
lm <- cbind(ids, lm)
# ordering measures by cluster labels
lm <- lm[order(lm$`labels$groups`), ]
lm$x <- NULL

perms <- data.frame(matrix(0, nrow = numboot, ncol = nsub))
bn <- nsub*bootsize #number of subjects in each bootstrap

##creating matrices to keep track of accuracies
acc_overall <- data.frame(matrix(0, nrow = 7, ncol = numboot))
acc_class <- data.frame(matrix(0, nrow = 1, ncol = 11))
list_adjrandindex <- data.frame(matrix(NA, nrow = numboot, ncol = 1)) 

### Running the actual random forest classifer
for(idx in 1:numboot){
  print(idx) # permutation number
  subjects <- t(train[idx, ]) # getting subjects for that permutation
  
  training_set <- cbind(subjects, lm)
  training_set <- training_set[which(training_set[,1] == "1"),] ## only getting the people for this perm
  ## shuffling labels
  training_set[,2] <- sample(training_set$`labels$groups`, replace = FALSE)
  labels_train <- as.factor(training_set[ ,2])
  training_set[,1:2] <- NULL
  
  ## setting up testing dataset - don't want to shuffle these
  subjects <- t(test[idx, ]) # getting subjects for that permutation
  testing_set <- cbind(subjects, lm)
  testing_set <- testing_set[which(testing_set[,1] == "1"),]
  labels_test <- as.factor(testing_set[ ,2])
  testing_set[,1:2] <- NULL
  # running the model
  model.rf = randomForest(labels_train ~. , data=training_set, ntree=501, mtry=5, importance=TRUE)
  
  rfpredict <- predict(model.rf, testing_set)
  stats <- confusionMatrix(rfpredict, labels_test) 
  
  acc_bygroup_perm <- as.data.frame(stats[4][1])
  names(acc_class) <- names(acc_bygroup_perm)
  acc_class <- rbind(acc_class, acc_bygroup_perm)
  
  row.names(acc_overall) <- row.names(as.data.frame(stats[3][1]))
  acc_overall[ , idx] <- as.data.frame(stats[3][1])
  list_adjrandindex[idx, 1] <- adj.rand.index(rfpredict, labels_test)
}

# calculating a significance cut-off (p<0.05) for adjusted rand index 
names(list_adjrandindex) <- c("index")
adjrandindex <- sort(list_adjrandindex$index)
quantile(adjrandindex, 0.95)
# calculating a significant (p<0.05) overall model accuracy
acc_overall <- as.data.frame(t(acc_overall))
acc_all <- acc_overall
acc_all <- sort(acc_all$Accuracy)
quantile(acc_all, 0.95)
# calculating a significance cut-off (p<0.05) for sensitivity and specificity for each group
acc_class <-acc_class[-1,] 
acc_class$label <- rep(1:4, times=100, each=1)
acc_class$label <- as.factor(acc_class$label)
acc_1 <- acc_class[which(acc_class$label == "1"), ]
acc_2 <- acc_class[which(acc_class$label == "2"), ]
acc_3 <- acc_class[which(acc_class$label == "3"), ]
acc_4 <- acc_class[which(acc_class$label == "4"), ]
# sensitivity
acc_1 <- acc_1[order(acc_1$byClass.Sensitivity), ]
acc_2 <- acc_2[order(acc_2$byClass.Sensitivity), ]
acc_3 <- acc_3[order(acc_3$byClass.Sensitivity), ]
acc_4 <- acc_4[order(acc_4$byClass.Sensitivity), ]
quantile(acc_1$byClass.Sensitivity, 0.95)
quantile(acc_2$byClass.Sensitivity, 0.95)
quantile(acc_3$byClass.Sensitivity, 0.95)
quantile(acc_4$byClass.Sensitivity, 0.95)
# specificity
acc_1 <- acc_class[which(acc_class$label == "1"), ]
acc_2 <- acc_class[which(acc_class$label == "2"), ]
acc_3 <- acc_class[which(acc_class$label == "3"), ]
acc_4 <- acc_class[which(acc_class$label == "4"), ]
acc_1 <- acc_1[order(acc_1$byClass.Specificity), ]
acc_2 <- acc_2[order(acc_2$byClass.Specificity), ]
acc_3 <- acc_3[order(acc_3$byClass.Specificity), ]
acc_4 <- acc_4[order(acc_4$byClass.Specificity), ]
quantile(acc_1$byClass.Specificity, 0.95)
quantile(acc_2$byClass.Specificity, 0.95)
quantile(acc_3$byClass.Specificity, 0.95)
quantile(acc_4$byClass.Specificity, 0.95)





