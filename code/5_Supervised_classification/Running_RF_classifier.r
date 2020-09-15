## Running a random forest classifer on POND SNF groups based on subsets of features
#1) All features
#2) Top 2 features from each data type
#3) Top 35 features
#4) Top 10 features

library(randomForest)
library(caret)
library(data.table)

directory <- ("output/5_Supervised_classification")

## inporting cluster labels
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

## gets the resampling subgroups 
source("code/5_Supervised_classification/setup_classifer_script.R")

numboot <- 100 # number of permutations 
nsub <- 176 #number of subjects
bootsize <- 0.8 #what percentage of participants to take per permuation

## getting the testing and training set
both <- dividing_data(numboot=numboot, nsub=nsub, bootsize=bootsize, labels=labels, measures=measures, ids=ids)
train <- both[1:numboot, ]
test <- both[(1+numboot):(2*numboot), ]

###########################################################################################
## Choosing which features will be included in the model

# All features included in the classification model

names <- read.csv("output/1_SNF_analysis/all_scores_k18_0.8_1000perms.csv", header = TRUE)
names <- names[1:135, 3]
#adjusting tract names
names <- gsub('-', '.', names)
names <- gsub('/', '.', names)
#setdiff(colnames(measures), names)
names <- as.data.frame(names)

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
importance_features <- data.frame(matrix(NA, nrow = length(measures), ncol = numboot))
importance_values <- data.frame(matrix(NA, nrow = length(measures), ncol = numboot))

### Running the actual random forest classifer
for(idx in 1:numboot){
  print(idx) # permutation number
  subjects <- t(train[idx, ]) # getting training set subjects for that permutation
  training_set <- cbind(subjects, lm) # adding cluster labels
  training_set <- training_set[which(training_set[,1] == "1"),] ## only getting the people for this perm
  labels_train <- as.factor(training_set[ ,2]) ## recording these labels
  training_set[,1:2] <- NULL
  
  ## setting up testing dataset
  subjects <- t(test[idx, ]) # getting testing subjects for that permutation
  testing_set <- cbind(subjects, lm)
  testing_set <- testing_set[which(testing_set[,1] == "1"),]
  labels_test <- as.factor(testing_set[ ,2])
  testing_set[,1:2] <- NULL
  
  ## training random forest
  model.rf = randomForest(labels_train ~. , data=training_set, ntree=501, mtry=80, importance=TRUE)

  #plot(svmfit, training_set)
  rfpredict <- predict(model.rf, testing_set)
  plot(rfpredict, labels_test)
  summary(rfpredict)
  
  # calculating the most important features
  import <- as.data.frame(importance(model.rf, type=1))
  import$Features <- row.names(import)
  import <- import[order(-import$MeanDecreaseAccuracy), ]
  importance_features[ ,idx] <- import$Features
  importance_values[ ,idx] <- import$MeanDecreaseAccuracy
  
  varImpPlot(model.rf, sort = TRUE)
  
  ### testing accuracy
  stats <- confusionMatrix(rfpredict, labels_test) 
  
  acc_bygroup_perm <- as.data.frame(stats[4][1])
  names(acc_class) <- names(acc_bygroup_perm)
  acc_class <- rbind(acc_class, acc_bygroup_perm)
  
  row.names(acc_overall) <- row.names(as.data.frame(stats[3][1]))
  acc_overall[ , idx] <- as.data.frame(stats[3][1])
  list_adjrandindex[idx, 1] <- adj.rand.index(rfpredict, labels_test)
}

model <- "all_mtry80"

write.csv(list_adjrandindex, file=file.path(directory, paste("list_adjrandidx_100perms_", model, ".csv", sep="")))
write.csv(acc_class, file=file.path(directory, paste("acc_by_group_100perms_", model, ".csv", sep="")))
write.csv(acc_overall, file=file.path(directory, paste("acc_overall_100perms_", model, ".csv", sep="")))

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
importance_features <- data.frame(matrix(NA, nrow = length(measures), ncol = numboot))
importance_values <- data.frame(matrix(NA, nrow = length(measures), ncol = numboot))

### Running the actual random forest classifer
for(idx in 1:numboot){
  print(idx) # permutation number
  subjects <- t(train[idx, ]) # getting training set subjects for that permutation
  training_set <- cbind(subjects, lm) # adding cluster labels
  training_set <- training_set[which(training_set[,1] == "1"),] ## only getting the people for this perm
  labels_train <- as.factor(training_set[ ,2]) ## recording these labels
  training_set[,1:2] <- NULL
  
  ## setting up testing dataset
  subjects <- t(test[idx, ]) # getting testing subjects for that permutation
  testing_set <- cbind(subjects, lm)
  testing_set <- testing_set[which(testing_set[,1] == "1"),]
  labels_test <- as.factor(testing_set[ ,2])
  testing_set[,1:2] <- NULL
  
  ## training random forest
  model.rf = randomForest(labels_train ~. , data=training_set, ntree=501, mtry=4, importance=TRUE)
  
  #plot(svmfit, training_set)
  rfpredict <- predict(model.rf, testing_set)
  plot(rfpredict, labels_test)
  summary(rfpredict)
  
  # calculating the most important features
  import <- as.data.frame(importance(model.rf, type=1))
  import$Features <- row.names(import)
  import <- import[order(-import$MeanDecreaseAccuracy), ]
  importance_features[ ,idx] <- import$Features
  importance_values[ ,idx] <- import$MeanDecreaseAccuracy
 
  varImpPlot(model.rf, sort = TRUE)
  
  ### testing accuracy
  stats <- confusionMatrix(rfpredict, labels_test) 
  
  acc_bygroup_perm <- as.data.frame(stats[4][1])
  names(acc_class) <- names(acc_bygroup_perm)
  acc_class <- rbind(acc_class, acc_bygroup_perm)
  
  row.names(acc_overall) <- row.names(as.data.frame(stats[3][1]))
  acc_overall[ , idx] <- as.data.frame(stats[3][1])
  list_adjrandindex[idx, 1] <- adj.rand.index(rfpredict, labels_test)
}

model <- "top2_mtry4"

write.csv(list_adjrandindex, file=file.path(directory, paste("list_adjrandidx_100perms_", model, ".csv", sep="")))
write.csv(acc_class, file=file.path(directory, paste("acc_by_group_100perms_", model, ".csv", sep="")))
write.csv(acc_overall, file=file.path(directory, paste("acc_overall_100perms_", model, ".csv", sep="")))

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
importance_features <- data.frame(matrix(NA, nrow = length(measures), ncol = numboot))
importance_values <- data.frame(matrix(NA, nrow = length(measures), ncol = numboot))

### Running the actual random forest classifer
for(idx in 1:numboot){
  print(idx) # permutation number
  subjects <- t(train[idx, ]) # getting training set subjects for that permutation
  training_set <- cbind(subjects, lm) # adding cluster labels
  training_set <- training_set[which(training_set[,1] == "1"),] ## only getting the people for this perm
  labels_train <- as.factor(training_set[ ,2]) ## recording these labels
  training_set[,1:2] <- NULL
  
  ## setting up testing dataset
  subjects <- t(test[idx, ]) # getting testing subjects for that permutation
  testing_set <- cbind(subjects, lm)
  testing_set <- testing_set[which(testing_set[,1] == "1"),]
  labels_test <- as.factor(testing_set[ ,2])
  testing_set[,1:2] <- NULL
  
  ## training random forest
  model.rf = randomForest(labels_train ~. , data=training_set, ntree=501, mtry=15, importance=TRUE)
  
  #plot(svmfit, training_set)
  rfpredict <- predict(model.rf, testing_set)
  plot(rfpredict, labels_test)
  summary(rfpredict)
  
  # calculating the most important features
  import <- as.data.frame(importance(model.rf, type=1))
  import$Features <- row.names(import)
  import <- import[order(-import$MeanDecreaseAccuracy), ]
  importance_features[ ,idx] <- import$Features
  importance_values[ ,idx] <- import$MeanDecreaseAccuracy
  
  varImpPlot(model.rf, sort = TRUE)
  
  ### testing accuracy
  stats <- confusionMatrix(rfpredict, labels_test) 
  
  acc_bygroup_perm <- as.data.frame(stats[4][1])
  names(acc_class) <- names(acc_bygroup_perm)
  acc_class <- rbind(acc_class, acc_bygroup_perm)
  
  row.names(acc_overall) <- row.names(as.data.frame(stats[3][1]))
  acc_overall[ , idx] <- as.data.frame(stats[3][1])
  list_adjrandindex[idx, 1] <- adj.rand.index(rfpredict, labels_test)
}

model <- "top35_15mtry"

write.csv(list_adjrandindex, file=file.path(directory, paste("list_adjrandidx_100perms_", model, ".csv", sep="")))
write.csv(acc_class, file=file.path(directory, paste("acc_by_group_100perms_", model, ".csv", sep="")))
write.csv(acc_overall, file=file.path(directory, paste("acc_overall_100perms_", model, ".csv", sep="")))


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
importance_features <- data.frame(matrix(NA, nrow = length(measures), ncol = numboot))
importance_values <- data.frame(matrix(NA, nrow = length(measures), ncol = numboot))

### Running the actual random forest classifer
for(idx in 1:numboot){
  print(idx) # permutation number
  subjects <- t(train[idx, ]) # getting training set subjects for that permutation
  training_set <- cbind(subjects, lm) # adding cluster labels
  training_set <- training_set[which(training_set[,1] == "1"),] ## only getting the people for this perm
  labels_train <- as.factor(training_set[ ,2]) ## recording these labels
  training_set[,1:2] <- NULL
  
  ## setting up testing dataset
  subjects <- t(test[idx, ]) # getting testing subjects for that permutation
  testing_set <- cbind(subjects, lm)
  testing_set <- testing_set[which(testing_set[,1] == "1"),]
  labels_test <- as.factor(testing_set[ ,2])
  testing_set[,1:2] <- NULL
  
  ## training random forest
  model.rf = randomForest(labels_train ~. , data=training_set, ntree=501, mtry=5, importance=TRUE)
  
  #plot(svmfit, training_set)
  rfpredict <- predict(model.rf, testing_set)
  plot(rfpredict, labels_test)
  summary(rfpredict)
  
  # calculating the most important features
  import <- as.data.frame(importance(model.rf, type=1))
  import$Features <- row.names(import)
  import <- import[order(-import$MeanDecreaseAccuracy), ]
  importance_features[ ,idx] <- import$Features
  importance_values[ ,idx] <- import$MeanDecreaseAccuracy
  
  varImpPlot(model.rf, sort = TRUE)
  
  ### testing accuracy
  stats <- confusionMatrix(rfpredict, labels_test) 
  
  acc_bygroup_perm <- as.data.frame(stats[4][1])
  names(acc_class) <- names(acc_bygroup_perm)
  acc_class <- rbind(acc_class, acc_bygroup_perm)
  
  row.names(acc_overall) <- row.names(as.data.frame(stats[3][1]))
  acc_overall[ , idx] <- as.data.frame(stats[3][1])
  list_adjrandindex[idx, 1] <- adj.rand.index(rfpredict, labels_test)
}

model <- "top10_mtry5"

write.csv(list_adjrandindex, file=file.path(directory, paste("list_adjrandidx_100perms_", model, ".csv", sep="")))
write.csv(acc_class, file=file.path(directory, paste("acc_by_group_100perms_", model, ".csv", sep="")))
write.csv(acc_overall, file=file.path(directory, paste("acc_overall_100perms_", model, ".csv", sep="")))









