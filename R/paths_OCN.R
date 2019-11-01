
paths_OCN <- function(OCN,
                      pathsRNlevel=FALSE){
  
  if (!("RN" %in% names(OCN))){
    stop('Missing fields in OCN. You should run aggregate_OCN prior to paths_OCN.')
  }
  
  # RN level
  RN_DownstreamPathLength <- spam(0,OCN$RN$Nnodes,OCN$RN$Nnodes)
  indices <- matrix(0,1000*OCN$RN$Nnodes,2)
  values <- numeric(1000*OCN$RN$Nnodes)
  RN_DownstreamPath <- vector("list",OCN$RN$Nnodes)
  counter <- 1
  for (i in 1:OCN$RN$Nnodes){RN_DownstreamPath[[i]] <- vector("list",OCN$RN$Nnodes)}
  for (i in 1:OCN$RN$Nnodes){
    Ups <- OCN$RN$Upstream[[i]]
    for (j in 1:length(Ups)){
      k <- Ups[j]
      Path <- k
      while (k != i) {
        k <- OCN$RN$DownNode[k]
        Path <- c(Path, k)
      }
      RN_DownstreamPath[[Ups[j]]][[i]] <- Path
      indices[counter, ] <- c(Ups[j],i)
      values[counter] <- sum(OCN$RN$Length[Path])
      counter <- counter + 1
      #RN_DownstreamPathLength[Ups[j],i] <- sum(OCN$RN$Length[Path])
    }
    # if ((i %% round(OCN$RN$Nnodes/100))==0){
    #   print(sprintf('paths at RN level: %0.1f%% completed',i*100/OCN$RN$Nnodes))
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
  
  # RN_DwnstrLength_unconnected <- matrix(data=0,nrow=OCN$RN$Nnodes,ncol=OCN$RN$Nnodes)
  # CountPaths <- rep(0,OCN$RN$Nnodes)
  # for (i in 1:OCN$RN$Nnodes){
  #   for (j in 1:OCN$RN$Nnodes){
  #     # check whether i and j belong to the same catchment
  #     if (OCN$RN$to_CM[i]==OCN$RN$to_CM[j]){
  #       k  <- intersect(RN_DownstreamPath[[i]][[OCN$RN$Outlet[OCN$RN$to_CM[i]]]],RN_DownstreamPath[[j]][[OCN$RN$Outlet[OCN$RN$to_CM[j]]]])[1]
  #       RN_DwnstrLength_unconnected[i,j] <- RN_DownstreamPathLength[i,k]
  #       #CountPaths[k] <- CountPaths[k]  + 1 
  #     }
  #   }
  # }
  # RN_DwnstrLength_unconnected <- RN_DwnstrLength_unconnected * (RN_DownstreamPathLength==0) # remove flow-connected reaches
  # RN_DwnstrLength_unconnected <- as(RN_DwnstrLength_unconnected,"sparseMatrix")
  # #print(sprintf('Elapsed time %.2f s',difftime(Sys.time(),t1,units='secs')),quote=FALSE)
  # #t1 <- Sys.time()
  
  # RN_DwnstrLength_unconnected <- matrix(data=0,nrow=OCN$RN$Nnodes,ncol=OCN$RN$Nnodes)
  # CountPaths <- rep(0,OCN$RN$Nnodes)
  # for (i in 1:OCN$RN$Nnodes){
  #   for (j in 1:OCN$RN$Nnodes){
  #     # check whether i and j belong to the same catchment
  #     if (OCN$RN$to_CM[i]==OCN$RN$to_CM[j]){
  #       k  <- intersect(RN_DownstreamPath[[i]][[OCN$RN$Outlet[OCN$RN$to_CM[i]]]],RN_DownstreamPath[[j]][[OCN$RN$Outlet[OCN$RN$to_CM[j]]]])[1]
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
  AG_DownstreamPathLength <- spam(0,OCN$AG$Nnodes,OCN$AG$Nnodes)
  AG_DownstreamPath <- vector("list",OCN$AG$Nnodes)
  indices <- matrix(0,1000*OCN$AG$Nnodes,2)
  values <- numeric(1000*OCN$AG$Nnodes)
  counter <- 1
  for (i in 1:OCN$AG$Nnodes){AG_DownstreamPath[[i]] <- vector("list",OCN$AG$Nnodes)}
  for (i in 1:OCN$AG$Nnodes){
    Ups <- OCN$AG$Upstream[[i]]
    for (j in 1:length(Ups)){
      k <- Ups[j]
      Path <- k
      while (k != i) {
        k <- OCN$AG$DownNode[k]
        Path <- c(Path, k)
      }
      AG_DownstreamPath[[Ups[j]]][[i]] <- Path
      #AG_DownstreamPathLength[Ups[j],i] <- sum(OCN$AG$Length[Path])
      indices[counter, ] <- c(Ups[j],i)
      values[counter] <- sum(OCN$AG$Length[Path])
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
  
  AG_DwnstrLength_unconnected <- spam(0,OCN$AG$Nnodes,OCN$AG$Nnodes)
  indices <- matrix(0,1000*OCN$AG$Nnodes,2)
  values <- numeric(1000*OCN$AG$Nnodes)
  counter <- 1
  for (i in 1:OCN$AG$Nnodes){
    for (j in 1:OCN$AG$Nnodes){
      # check whether i and j belong to the same catchment
      if (OCN$AG$to_CM[i]==OCN$AG$to_CM[j] && AG_DownstreamPathLength[i,j]==0){
        k  <- intersect(AG_DownstreamPath[[i]][[OCN$AG$Outlet[OCN$AG$to_CM[i]]]],AG_DownstreamPath[[j]][[OCN$AG$Outlet[OCN$AG$to_CM[j]]]])[1]
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
  if (pathsRNlevel==TRUE){ OCN$RN$DownstreamPath <- RN_DownstreamPath}
  OCN$RN$DownstreamPathLength <- RN_DownstreamPathLength
  # OCN$RN$DwnstrLength_unconnected <- RN_DwnstrLength_unconnected
  OCN$AG$DownstreamPath <- AG_DownstreamPath
  OCN$AG$DownstreamPathLength <- AG_DownstreamPathLength
  OCN$AG$DwnstrLength_unconnected <- AG_DwnstrLength_unconnected
  
  
  return(OCN)
  
}


