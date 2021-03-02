
paths_OCN <- function(OCN,
                      pathsRN=FALSE,
                      includeDownstreamNode=FALSE){
  
  if (!("RN" %in% names(OCN))){
    stop('Missing fields in OCN. You should run aggregate_OCN prior to paths_OCN.')
  }
  
  # RN level
  RN_DownstreamPathLength <- spam(0,OCN$RN$nNodes,OCN$RN$nNodes)
  indices <- matrix(0,OCN$RN$nNodes^2,2)
  values <- numeric(OCN$RN$nNodes^2)
  RN_DownstreamPath <- vector("list",OCN$RN$nNodes)
  counter <- 1
  for (i in 1:OCN$RN$nNodes){RN_DownstreamPath[[i]] <- vector("list",OCN$RN$nNodes)}
  for (i in 1:OCN$RN$nNodes){
    Ups <- OCN$RN$upstream[[i]]
    for (j in 1:length(Ups)){
      k <- Ups[j]
      Path <- k
      while (k != i) {
        k <- OCN$RN$downNode[k]
        Path <- c(Path, k)
      }
      RN_DownstreamPath[[Ups[j]]][[i]] <- Path
      indices[counter, ] <- c(Ups[j],i)
      values[counter] <- sum(OCN$RN$leng[Path]) 
      if (includeDownstreamNode==FALSE){
        values[counter] <- values[counter] - OCN$RN$leng[k]
      }
      counter <- counter + 1
      #RN_DownstreamPathLength[Ups[j],i] <- sum(OCN$RN$leng[Path])
    }
    # if ((i %% round(OCN$RN$nNodes/100))==0){
    #   print(sprintf('paths at RN level: %0.1f%% completed',i*100/OCN$RN$nNodes))
    # }
  }
  indices <- indices[1:(counter-1),]
  values <- values[1:(counter-1)]
  RN_DownstreamPathLength[indices] <- values
  
  #RN_DownstreamPathLength <- as(RN_DownstreamPathLength,"sparseMatrix")
  #print(sprintf('Elapsed time %.2f s',difftime(Sys.time(),t1,units='secs')),quote=FALSE) 
  #t1 <- Sys.time()
  
  #print('Downstream/Upstream path lengths for unconnected reaches...',quote=FALSE) 
  # DwnstrLength_unconnected[i,j] is the length of the downstream path joining i to k, where k is the junction node between two flow-unconnected reaches i and j 
  # CountPaths[k] is the number of paths connecting pairs of reaches passing through k, considering also k as a start/end point (i.e. the path connecting k to k is counted)
  
  # RN_DwnstrLength_unconnected <- matrix(data=0,nrow=OCN$RN$nNodes,ncol=OCN$RN$nNodes)
  # CountPaths <- rep(0,OCN$RN$nNodes)
  # for (i in 1:OCN$RN$nNodes){
  #   for (j in 1:OCN$RN$nNodes){
  #     # check whether i and j belong to the same catchment
  #     if (OCN$RN$toCM[i]==OCN$RN$toCM[j]){
  #       k  <- intersect(RN_DownstreamPath[[i]][[OCN$RN$outlet[OCN$RN$toCM[i]]]],RN_DownstreamPath[[j]][[OCN$RN$outlet[OCN$RN$toCM[j]]]])[1]
  #       RN_DwnstrLength_unconnected[i,j] <- RN_DownstreamPathLength[i,k]
  #       #CountPaths[k] <- CountPaths[k]  + 1 
  #     }
  #   }
  # }
  # RN_DwnstrLength_unconnected <- RN_DwnstrLength_unconnected * (RN_DownstreamPathLength==0) # remove flow-connected reaches
  # RN_DwnstrLength_unconnected <- as(RN_DwnstrLength_unconnected,"sparseMatrix")
  # #print(sprintf('Elapsed time %.2f s',difftime(Sys.time(),t1,units='secs')),quote=FALSE)
  # #t1 <- Sys.time()
  
  # RN_DwnstrLength_unconnected <- matrix(data=0,nrow=OCN$RN$nNodes,ncol=OCN$RN$nNodes)
  # CountPaths <- rep(0,OCN$RN$nNodes)
  # for (i in 1:OCN$RN$nNodes){
  #   for (j in 1:OCN$RN$nNodes){
  #     # check whether i and j belong to the same catchment
  #     if (OCN$RN$toCM[i]==OCN$RN$toCM[j]){
  #       k  <- intersect(RN_DownstreamPath[[i]][[OCN$RN$outlet[OCN$RN$toCM[i]]]],RN_DownstreamPath[[j]][[OCN$RN$outlet[OCN$RN$toCM[j]]]])[1]
  #       RN_DwnstrLength_unconnected[i,j] <- RN_DownstreamPathLength[i,k]
  #       #CountPaths[k] <- CountPaths[k]  + 1 
  #     }
  #   }
  # }
  # RN_DwnstrLength_unconnected <- RN_DwnstrLength_unconnected * (RN_DownstreamPathLength==0) # remove flow-connected reaches
  # RN_DwnstrLength_unconnected <- as(RN_DwnstrLength_unconnected,"sparseMatrix")
  # #print(sprintf('Elapsed time %.2f s',difftime(Sys.time(),t1,units='secs')),quote=FALSE)
  # #t1 <- Sys.time()
  
  # AG level
  AG_DownstreamPathLength <- spam(0,OCN$AG$nNodes,OCN$AG$nNodes)
  AG_DownstreamPath <- vector("list",OCN$AG$nNodes)
  indices <- matrix(0,OCN$AG$nNodes^2,2)
  values <- numeric(OCN$AG$nNodes^2)
  counter <- 1
  for (i in 1:OCN$AG$nNodes){AG_DownstreamPath[[i]] <- vector("list",OCN$AG$nNodes)}
  for (i in 1:OCN$AG$nNodes){
    Ups <- OCN$AG$upstream[[i]]
    for (j in 1:length(Ups)){
      k <- Ups[j]
      Path <- k
      while (k != i) {
        k <- OCN$AG$downNode[k]
        Path <- c(Path, k)
      }
      AG_DownstreamPath[[Ups[j]]][[i]] <- Path
      #AG_DownstreamPathLength[Ups[j],i] <- sum(OCN$AG$leng[Path])
      indices[counter, ] <- c(Ups[j],i)
      values[counter] <- sum(OCN$AG$leng[Path]) 
      if (includeDownstreamNode==FALSE){
        values[counter] <- values[counter]  - OCN$AG$leng[k]
      }
      counter <- counter + 1
    }
  }
  indices <- indices[1:(counter-1), ]
  values <- values[1:(counter-1)]
  AG_DownstreamPathLength[indices] <- values
  #AG_DownstreamPathLength <- as(AG_DownstreamPathLength,"sparseMatrix")
  #print(sprintf('Elapsed time %.2f s',difftime(Sys.time(),t1,units='secs')),quote=FALSE) 
  #t1 <- Sys.time()
  
  #print('Downstream/Upstream path lengths for unconnected reaches...',quote=FALSE) 
  # DwnstrLength_unconnected[i,j] is the length of the downstream path joining i to k, where k is the junction node between two flow-unconnected reaches i and j 
  # CountPaths[k] is the number of paths connecting pairs of reaches passing through k, considering also k as a start/end point (i.e. the path connecting k to k is counted)
  
  AG_DwnstrLength_unconnected <- spam(0,OCN$AG$nNodes,OCN$AG$nNodes)
  indices <- matrix(0,OCN$AG$nNodes^2,2)
  values <- numeric(OCN$AG$nNodes^2)
  counter <- 1
  for (i in 1:OCN$AG$nNodes){
    for (j in 1:OCN$AG$nNodes){
      # check whether i and j belong to the same catchment
      if (OCN$AG$toCM[i]==OCN$AG$toCM[j] && AG_DownstreamPathLength[i,j]==0){
        k  <- intersect(AG_DownstreamPath[[i]][[OCN$AG$outlet[OCN$AG$toCM[i]]]],AG_DownstreamPath[[j]][[OCN$AG$outlet[OCN$AG$toCM[j]]]])[1]
        #AG_DwnstrLength_unconnected[i,j] <- AG_DownstreamPathLength[i,k]
        indices[counter,] <- c(i,j)
        values[counter] <- AG_DownstreamPathLength[i,k]
        counter <- counter + 1
        #CountPaths[k] <- CountPaths[k]  + 1 
      }
    }
  }
  indices <- indices[1:(counter-1), ]
  values <- values[1:(counter-1)]
  AG_DwnstrLength_unconnected[indices] <- values
  
  #AG_DwnstrLength_unconnected <- AG_DwnstrLength_unconnected * (AG_DownstreamPathLength==0) # remove flow-connected reaches
  #AG_DwnstrLength_unconnected <- as(AG_DwnstrLength_unconnected,"sparseMatrix")
  #print(sprintf('Elapsed time %.2f s',difftime(Sys.time(),t1,units='secs')),quote=FALSE)
  #t1 <- Sys.time()
  
  ## copy variables into OCN
  if (pathsRN==TRUE){ OCN$RN$downstreamPath <- RN_DownstreamPath}
  OCN$RN$downstreamPathLength <- RN_DownstreamPathLength
  # OCN$RN$DwnstrLength_unconnected <- RN_DwnstrLength_unconnected
  OCN$AG$downstreamPath <- AG_DownstreamPath
  OCN$AG$downstreamPathLength <- AG_DownstreamPathLength
  OCN$AG$downstreamLengthUnconnected <- AG_DwnstrLength_unconnected
  
  
  invisible(OCN)
  
}


