###################
#     Robust Core-clustering function
###################

# This set of functions finds the share of times a random draw of your 
# sample will put each pair of individuals in the same cluster. 
# This allows you to (1) only include individuals who reliably cluster within
# their own cluster in the analysis and (2) identify cluster that stably 
# persist given a certain number of samplings of your data. 

## Finds the core clustering matrices. 
## (1) The first matrix it finds is nxn with each cell representing the share of times 
## each individual i clusters with individual j. 
## (2) The second matrix sets values that represent below random shares of cluster 
## (i.e. share < (1/number of clusters)^2) to zero

RobustCoreClusteringMatrix <- function(feature.affinity.mat.list,num.clusts,exp.num.samples = 1000,percent.sample = 0.8, clust.type = 2,seed = 3333){
  ## feature.affinity.mat.list = List of features you used to generate your SNF matrix 
    ## (same as input to the SNF function)
  
  ## num.clusts = the number of clusters you'd like to test the co-clustering of
  ## exp.num.samples = the expected number of times you would like to sample each individual
  ## percent.sample = the percentage of individuals you'd like to sample each time
  ## clust.type = spectral clustering type 1 or 2 (default is 2, should likely never change)
  ## seed = a seed so every time you run the loop, it's the same (same random draw)
  
  n.patients = nrow(feature.affinity.mat.list[[1]])  
  
  for(m in 1:length(feature.affinity.mat.list)){
    rownames(feature.affinity.mat.list[[m]]) = colnames(feature.affinity.mat.list[[m]]) = as.character(1:n.patients)  
  }
  
  ## Creating empty matrix to put number of times the individuals are *sampled* together
  n.pair.clust.mat <- matrix(0,nrow=n.patients,ncol=n.patients)
  ## Naming the rows and columns '1','2',...
  colnames(n.pair.clust.mat) <- rownames(n.pair.clust.mat) <- as.character(1:n.patients)
  
  ## Create an empty matrix to count the number of times a pair of individuals *cluster* together
  count.coclust.mat <- matrix(0,nrow=n.patients,ncol=n.patients)
  ## Naming the rows and columns '1','2',...
  colnames(count.coclust.mat) <- rownames(count.coclust.mat) <- as.character(1:n.patients)

  set.seed(seed)

  if(clust.type == 2){
    for(k in 1:(exp.num.samples*percent.sample)){
      ## Sample <percentage> of your patients
      sample.set <- sample(x = 1:n.patients,size = round(percent.sample*n.patients))
      
      ## Subset those patients from your separate data types
      sample.dat.mat.list <- lapply(feature.affinity.mat.list,function(x)x[as.character(sample.set),as.character(sample.set)])
      
      ## Run SNF on your sub sample of individuals from above
      sample.snf <- SNF(Wall = sample.dat.mat.list,K = K,t = T)
      ## Name rows and columns
      colnames(sample.snf) <- rownames(sample.snf) <- as.character(sample.set)
      #   str(sample.snf)
      
      n.clust <- num.clusts
      
      ## Run spectral clustering on SNF matrix using clusters from above
      sample.groups <- spectralClustering2(W = sample.snf,C = n.clust)      

      ## Create a dataframe making sub sampled id's correspond to their group assignment
      sample.group.df <- data.frame(cbind(as.character(sample.set),sample.groups))
      
      ## Naming columns in dataframe created above
      names(sample.group.df) <- c('id','group.id')
      sample.group.df$id <- as.numeric(as.character(sample.group.df$id))
      #   head(sample.group.df)
      
      ## Adding 1 to each cell of matrix that corresponds to 2 people that co-cluster
      for(i in sample.set){
        ## outerloop for each individual sampled
        for(j in sample.set[sample.set != i]){       ## innerloop for each individual except for i (so we can compare the others all to i)
          #       n.pair.clust.mat[i,j] <- n.pair.clust.mat[j,i] <- n.pair.clust.mat[i,j] + 1 ## add 1 b/c they are both sampled
          n.pair.clust.mat[i,j] <- n.pair.clust.mat[i,j] + 1 ## add 1 b/c they are both sampled
          if(sample.group.df$group.id[sample.group.df$id == i] == sample.group.df$group.id[sample.group.df$id == j]){ ## if they co-cluster
            #         count.coclust.mat[i,j] <- count.coclust.mat[j,i] <- count.coclust.mat[i,j] + 1 ## we add one to this cell of the matrix corresponding to the pair
            count.coclust.mat[i,j] <- count.coclust.mat[i,j] + 1 ## we add one to this cell of the matrix corresponding to the pair
          }
        }
      }
      
      ## print share of iterations complete so you can keep track
      if(k%%((exp.num.samples/10)*percent.sample) == 0){    
        cat(100*(k/(exp.num.samples*percent.sample)),'% of iterations complete.\n')
      }
    }
    
  }

  else{
    for(k in 1:(exp.num.samples*percent.sample)){
      ## Sample <percentage> of your patients
      sample.set <- sample(x = 1:n.patients,size = round(percent.sample*n.patients))
      
      ## Subset those patients from your separate data types
      sample.dat.mat.list <- lapply(feature.affinity.mat.list,function(x)x[as.character(sample.set),as.character(sample.set)])
      
      ## Run SNF on your sub sample of individuals from above
      sample.snf <- SNF(Wall = sample.dat.mat.list,K = K,t = T)
      ## Name rows and columns
      colnames(sample.snf) <- rownames(sample.snf) <- as.character(sample.set)
      #   str(sample.snf)
      
      n.clust <- num.clusts
      
      ## Run spectral clustering on SNF matrix using clusters from above
      sample.groups <- spectralClustering(affinity = sample.snf,K = n.clust)
      
      ## Create a dataframe making sub sampled id's correspond to their group assignment
      sample.group.df <- data.frame(cbind(as.character(sample.set),sample.groups))
      
      ## Naming columns in dataframe created above
      names(sample.group.df) <- c('id','group.id')
      sample.group.df$id <- as.numeric(as.character(sample.group.df$id))
      #   head(sample.group.df)
      
      ## Adding 1 to each cell of matrix that corresponds to 2 people that co-cluster
      for(i in sample.set){
        ## outerloop for each individual sampled
        for(j in sample.set[sample.set != i]){       ## innerloop for each individual except for i (so we can compare the others all to i)
          #       n.pair.clust.mat[i,j] <- n.pair.clust.mat[j,i] <- n.pair.clust.mat[i,j] + 1 ## add 1 b/c they are both sampled
          n.pair.clust.mat[i,j] <- n.pair.clust.mat[i,j] + 1 ## add 1 b/c they are both sampled
          if(sample.group.df$group.id[sample.group.df$id == i] == sample.group.df$group.id[sample.group.df$id == j]){ ## if they co-cluster
            #         count.coclust.mat[i,j] <- count.coclust.mat[j,i] <- count.coclust.mat[i,j] + 1 ## we add one to this cell of the matrix corresponding to the pair
            count.coclust.mat[i,j] <- count.coclust.mat[i,j] + 1 ## we add one to this cell of the matrix corresponding to the pair
          }
        }
      }
      
      ## print share of iterations complete so you can keep track
      if(k%%((exp.num.samples/10)*percent.sample) == 0){    
        cat(100*(k/(exp.num.samples*percent.sample)),'% of iterations complete.\n')
      }
    }
    
  }
  ## Running co-clustering loop
  
  ### Divide matrix that counts co-clustering by matrix that counts number of times the pair were sampled together
  ## create empty matrix to put new values into
  share.coclust.mat <- matrix(ncol=n.patients,nrow=n.patients)
  ## name rows and columns of the matrix
  colnames(share.coclust.mat) <- rownames(share.coclust.mat) <- as.character(1:n.patients)
  ## iterate over every patient (outerloop) relative to every other patient (innerloop)
  for(i in as.character(1:n.patients)){
    for(j in as.character(1:n.patients)){
      ## Divide matrix element that counts co-clustering by matrix that counts number of times the pair were sampled together 
      share.coclust.mat[i,j] <- share.coclust.mat[j,i] <- count.coclust.mat[i,j]/n.pair.clust.mat[i,j]
    }
  }
  
  ## putting the average value in the diagonal in the matrix
  diag(share.coclust.mat) <- mean(na.omit(c(share.coclust.mat)))

  ### new.mat is the co-clustering mat with non-reliable/not above random
  ###   co-clustering pair values set to zero
  new.mat <- share.coclust.mat
  new.mat[which(share.coclust.mat < (1/num.clusts)^2)] <- 0

  return.list <- list(share.coclust.mat, new.mat)
  names(return.list) <- c('Dense Core Cluster Matrix','Sparse Core Cluster Matrix')

  return(return.list)
}

## This function returns a dataframe that gives the individual's id's and their group membership
## with an updated spectral clustering function
RobustCoreClusteringClusters <- function(core.clustering.list,num.clusts,id.vec = 1:nrow(core.clustering.list[[1]]), clust.type = 2, verbose = TRUE){
  require('SNFtool')

  ## core.clustering.list = the output from the RobustCoreClusteringMatrix function
    ## dense.core.cluster.matrix = 'Dense Core Cluster Matrix' from RobustCoreClusteringMatrix function
    ## sparse.core.cluster.matrix = 'Sparse Core Cluster Matrix' from RobustCoreClusteringMatrix function
  
  ## num.clusts = the number of clusters you'd like to test the co-clustering of
  ## id.vec = id's on the rows/columns of your affinity matrix
  ## clust.type = keep equal to 2 for consistent spectral clustering, any other value will 
    ## run standard/less consistent spectral clustering
  
  
  share.coclust.mat = core.clustering.list[[2]]
  
  n.inds <- nrow(share.coclust.mat)

  if(clust.type == 2){
    co.clust.groups <- spectralClustering2(W = share.coclust.mat,C = num.clusts)
    id.grp.df <- data.frame(cbind(as.character(id.vec),co.clust.groups))
    names(id.grp.df) <- c('id','groups')
    id.grp.df$id <- as.character(id.grp.df$id)
    
    return(id.grp.df)
    
    if(verbose == TRUE){
      for(m in 1:length(unique(id.grp.df$groups))){
        co.clust.sub <- share.coclust.mat[id.grp.df$id[id.grp.df$groups == m],id.grp.df$id[id.grp.df$groups == m]]
        n = nrow(id.grp.df[id.grp.df$groups == m,])
        n.reliable = (length(c(co.clust.sub)[c(co.clust.sub) > (1/num.clusts)^2])/2)
        n.pairs = (((n^2)-n)/2)
        share.reliable <- n.reliable/n.pairs
        cat('group',m,'\n')
        cat('size of group:',n,'\n')
        cat('number of non-random pairs',n.reliable,'\n')
        cat('total number of pairs',n.pairs,'\n')
        cat('share reliable above random chance',share.reliable,'\n')
      }    
    }      
  }
  else{
    co.clust.groups <- spectralClustering(affinity = share.coclust.mat,K = num.clusts)
    
    id.grp.df <- data.frame(cbind(as.character(id.vec),co.clust.groups))
    names(id.grp.df) <- c('id','groups')
    id.grp.df$id <- as.character(id.grp.df$id)
    
    return(id.grp.df)
    
    if(verbose == TRUE){
      for(m in 1:length(unique(id.grp.df$groups))){
        co.clust.sub <- share.coclust.mat[id.grp.df$id[id.grp.df$groups == m],id.grp.df$id[id.grp.df$groups == m]]
        n = nrow(id.grp.df[id.grp.df$groups == m,])
        n.reliable = (length(c(co.clust.sub)[c(co.clust.sub) > (1/num.clusts)^2])/2)
        n.pairs = (((n^2)-n)/2)
        share.reliable <- n.reliable/n.pairs
        cat('group',m,'\n')
        cat('size of group:',n,'\n')
        cat('number of non-random pairs',n.reliable,'\n')
        cat('total number of pairs',n.pairs,'\n')
        cat('share reliable above random chance',share.reliable,'\n')
      }    
    }  
  }
}

## The updated spectral clustering function
spectralClustering2=function(W,C){
  ## W = affinity matrix
  ## C = number of clusters
  
  require(stats)
  D=diag(apply(W,1,sum));U=D-W;L=U;
  "%^%" <- function(M,power)with(eigen(M), vectors %*% (values^power * solve(vectors)));
  L=(D%^%(-0.5))%*%U%*%(D%^%(-0.5))
  #L=diag(nrow(U))-solve(D)%*%W
  #L=solve(D)%*%U
  evl=eigen(L,symmetric=TRUE)
  ev=evl$vectors[,(ncol(evl$vectors)-C+1):ncol(evl$vectors)]
  lab=kmeans(ev,centers=C,iter.max=500,nstart=500)$cluster;
  #lab=apply(ev,1,which.max)
  #x=hclust(dist(ev)); lab=cutree(x,C) 
  return(lab)
}

## This functions returns individuals who do not reliably cluster with anyone in the cohort
## (this ideally should be zero)
UnrelClustPatientsFullCohort <- function(core.clustering.list,verbose = TRUE){
  ## core.clustering.list = list(dense.core.cluster.matrix,sparse.core.cluster.matrix)
    ## This is the output from the RobustCoreClusteringMatrix function
  
  sparse.mat <- core.clustering.list[[2]]
  n.inds <- ncol(sparse.mat)
  
  unreliable.inds.full.cohort <- c()

  k=1
  for(i in as.character(1:n.inds)){
    if(sum(sparse.mat[i,]) == mean(diag(sparse.mat))){
      unreliable.inds.full.cohort[k] <- i
      k = k + 1
    }
  }
  
  if(verbose == TRUE){
    cat(length(unreliable.inds.full.cohort),'individuals don\'t reliably cluster with anyone.')
  }
  return(unreliable.inds.full.cohort)  
}

## This function returns the individuals that don't reliably cluster with anyone within their own group 
## (this ideally should be zero)
UnrelClustPatientsByGrp <- function(core.clustering.list,id.group.df,verbose = TRUE){
  ## core.clustering.list = list(dense.core.cluster.matrix,sparse.core.cluster.matrix)
    ## This is the output from the RobustCoreClusteringMatrix function
  
  ## id.group.df = output from RobustCoreClusteringClusters function
    ## this should be a dataframe with id's in the first column, groups in the second column
    ## first column name = 'id', second column name = 'group'
  
  sparse.mat <- core.clustering.list[[2]]
  n.inds <- ncol(sparse.mat)
  
  unreliable.inds.by.group <- c()
  group <- c()
  k=1
  for(i in id.group.df$id){
    group.for.i <- id.group.df$group[as.character(id.group.df$id) == i]
    in.ind.i.group.not.ind.i <- (id.group.df$id[id.group.df$group == group.for.i])[(id.group.df$id[id.group.df$group == group.for.i]) != i]
    if(sum(c(sparse.mat[i,in.ind.i.group.not.ind.i])) == 0){
      unreliable.inds.by.group[k] <- i    
      group[k] <- group.for.i   
      k = k + 1
    }
  }

  if(length(unreliable.inds.by.group) == 0){
    cat('No individuals cluster unreliably within their given group.')
  }
  else{
    unrel.inds.grp.df <- data.frame(cbind(unreliable.inds.by.group,group))
    names(unrel.inds.grp.df) <- c("id","groups")
    
    return(unrel.inds.grp.df) ## these people are not a member of reliable pair within their group
  }
}
