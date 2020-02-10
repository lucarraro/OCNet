
find_area_threshold_OCN <- function(OCN,
                                    thrValues=seq(OCN$cellsize^2,min(OCN$CM$A),OCN$cellsize^2),
                                    maxReachLength=Inf,
                                    streamOrderType="Strahler",
                                    displayUpdates=0){
  
  if (!("slope" %in% names(OCN$FD))){
    stop('Missing fields in OCN. You should run landscape_OCN prior to find_area_threshold_OCN.')
  }
  
  if (max(thrValues) > min(OCN$CM$A)){
    stop('max(thrValues) cannot be larger than min(OCN$CM$A)')
  }
  if (length(thrValues)>5000 && OCN$FD$nNodes>5000){
    message('Running...\n', appendLF = FALSE)
    message('To reduce computational time, reduce the length of thrValues.\n', appendLF = FALSE)
  }
  
  vec_nNodesRN <- numeric(length(thrValues))
  vec_nNodesAG <- numeric(length(thrValues))
  vec_DrainageDensity <- numeric(length(thrValues))
  vec_StreamOrder <- numeric(length(thrValues))
  
  for (thr in 1:(length(thrValues))){
    
    if (thr < length(thrValues)){
    if (displayUpdates==1) {message(sprintf('\r%.2f%% completed',100*thr/length(thrValues)), appendLF = FALSE)}
    } else {
      if (displayUpdates==1) {message(sprintf('\r%.2f%% completed',100*thr/length(thrValues)))}  
    }
    # print(sprintf('Evaluating A_thr = %.0f m2',thrValues[thr]))
    
    RN_mask <- as.vector(OCN$FD$A >= thrValues[thr]) 
    FD_to_RN <- RN_mask*cumsum(as.numeric(RN_mask)) # FD_to_RN[i] is the pixel ID at the RN level of the pixel whose ID at the FD level is i
    # if pixel i at FD level doesn't belong to RN, then FD_to_RN[i]=0
    
    nNodes_RN <- sum(RN_mask)
    vec_nNodesRN[thr] <- nNodes_RN
    
    if (nNodes_RN > 1){
      
      W_RN <- OCN$FD$W[RN_mask,,drop=FALSE]
      W_RN <- W_RN[,RN_mask,drop=FALSE]
      
      Outlet_RN <- FD_to_RN[OCN$FD$outlet]
      Outlet_RN <- Outlet_RN[Outlet_RN!=0] # remove outlets if the corresponding catchment size is lower than A_threshold
      DownNode_RN <- numeric(nNodes_RN)
      tmp <- W_RN@rowpointers
      NotOutlet <- which((tmp[-1] - tmp[-length(tmp)])==1)
      DownNode_RN[NotOutlet] <- W_RN@colindices
      
      Length_RN <- OCN$FD$leng[RN_mask]
      
      # Drainage density
      vec_DrainageDensity[thr] <- sum(Length_RN)/(OCN$dimX*OCN$dimY*OCN$cellsize^2)
      
      # Connectivity indices at pixel level
      DegreeIn <- colSums(W_RN)
      DegreeOut <- rowSums(W_RN)
      Confluence <- DegreeIn>1
      Source <- DegreeIn==0
      SourceOrConfluence <- Source|Confluence
      ConfluenceNotOutlet <- Confluence&(DownNode_RN!=0)
      ChannelHeads <- SourceOrConfluence  #Source|ConfluenceNotOutlet
      
      OutletNotChannelHead <- (DownNode_RN==0)&(!ChannelHeads)
      IsNodeAG <- SourceOrConfluence|OutletNotChannelHead
      whichNodeAG <- which(IsNodeAG)
      
      
      nNodes_AG <- sum(IsNodeAG)
      Length_AG <- numeric(nNodes_AG)
      RN_to_AG <- numeric(nNodes_RN)
      reachID <- 1
      while (length(whichNodeAG) != 0){ # explore all AG Nodes
        i <- whichNodeAG[1] # select the first
        RN_to_AG[i] <- reachID; j <- DownNode_RN[i] 
        Length_AG[reachID] <- Length_RN[i]
        while (!IsNodeAG[j] && j!=0 && Length_AG[reachID] <= maxReachLength) {
          RN_to_AG[j] <- reachID 
          Length_AG[reachID] <-  Length_AG[reachID] + Length_RN[j]
          j_old <- j
          j <- DownNode_RN[j]} 
        if (Length_AG[reachID] > maxReachLength){
          j <- j_old
          Length_AG[reachID] <-  Length_AG[reachID] - Length_RN[j]
          ChannelHeads[j] <- 1
          whichNodeAG <- c(whichNodeAG,j)}
        reachID <- reachID + 1
        whichNodeAG <- whichNodeAG[-1]
      }
      nNodes_AG <- length(Length_AG)
      
      vec_nNodesAG[thr] <- nNodes_AG
      
      # RN_to_AG <- numeric(nNodes_RN)
      # reachID <- 1
      # for (i in 1:nNodes_RN){
      #   if (Source_Or_ConfluenceNotOutlet[i]) {
      #     RN_to_AG[i] <- reachID; j <- DownNode_RN[i]
      #      while (!Source_Or_ConfluenceNotOutlet[j] && j!=0) {
      #       RN_to_AG[j] <- reachID; j <- DownNode_RN[j]} 
      #     reachID <- reachID + 1
      #   }}
      
      #print('W matrix at AG level...',quote=FALSE); 
      # Adjacency matrix at reach level
      DownNode_AG <- numeric(vec_nNodesAG[thr])
      #W_AG <- sparseMatrix(i=1,j=1,x=0,dims=c(vec_nNodes[thr],vec_nNodes[thr]))
      W_AG <- spam(0,vec_nNodesAG[thr],vec_nNodesAG[thr])
      indices <- matrix(0,vec_nNodesAG[thr],2)
      
      for (i in 1:nNodes_RN){ 
        if (DownNode_RN[i] != 0 && RN_to_AG[DownNode_RN[i]] != RN_to_AG[i]) {
          DownNode_AG[RN_to_AG[i]] <- RN_to_AG[DownNode_RN[i]]
          #W_AG[RN_to_AG[i],DownNode_AG[RN_to_AG[i]]] <- 1
          indices[RN_to_AG[i],] <- c(RN_to_AG[i],DownNode_AG[RN_to_AG[i]])
        }
      }
      
      
      # for (i in 1:nNodes_RN){ 
      #   if (DownNode_RN[i]==0) {
      #     DownNode_AG[RN_to_AG[i]] <- reachID
      #     #W_AG[RN_to_AG[i],reachID] <- 1
      #     indices[RN_to_AG[i],] <- c(RN_to_AG[i],reachID)
      #     reachID <- reachID + 1
      #   } else if (RN_to_AG[DownNode_RN[i]] != RN_to_AG[i]) {
      #     DownNode_AG[RN_to_AG[i]] <- RN_to_AG[DownNode_RN[i]]
      #     #W_AG[RN_to_AG[i],DownNode_AG[RN_to_AG[i]]] <- 1
      #     indices[RN_to_AG[i],] <- c(RN_to_AG[i],DownNode_AG[RN_to_AG[i]])
      #   }}
      indices <- indices[-which(indices[,1]==0),]
      W_AG[indices] <- 1
      Outlet_AG <- which(DownNode_AG==0)
      
      # Upstream_AG : list containing IDs of all reaches upstream of each reach (plus reach itself)
      Upstream_AG <- vector("list",vec_nNodesAG[thr])
      Nupstream_AG <- numeric(vec_nNodesAG[thr])
      for (i in 1:vec_nNodesAG[thr]){
        UpOneLevel <- which(DownNode_AG==i) # find reaches at one level upstream
        Upstream_AG[[i]] <- UpOneLevel      # add them to the list
        while (length(UpOneLevel)!=0) { # continue until there are no more reaches upstream
          ContinuePath <- UpOneLevel # jump 1 level above
          UpOneLevel <- which(DownNode_AG %in% ContinuePath) # find reaches at one level upstream
          Upstream_AG[[i]] <- c(Upstream_AG[[i]],UpOneLevel) # add them to the list
        }
        Upstream_AG[[i]] <- c(Upstream_AG[[i]],i)
        Nupstream_AG[i] <- length(Upstream_AG[[i]])
      }
      
      
      ## Calculate stream Order
      if (streamOrderType=="Strahler"){
        # calculate Strahler stream order
        StreamOrder_AG <- numeric(vec_nNodesAG[thr])
        for (i in 1:vec_nNodesAG[thr]){
          j <- order(Nupstream_AG)[i] # index that explores reaches in a downstream direction
          tmp <- which(DownNode_AG==j) # set of reaches draining into j
          if (length(tmp)>0){
            IncreaseOrder <- sum(StreamOrder_AG[tmp]==max(StreamOrder_AG[tmp])) # check whether tmp reaches have the same stream order
            if (IncreaseOrder > 1) {
              StreamOrder_AG[j] <- 1 + max(StreamOrder_AG[tmp]) # if so, increase stream order
            } else {StreamOrder_AG[j] <- max(StreamOrder_AG[tmp])} # otherwise, keep previous stream order
          } else {StreamOrder_AG[j] <- 1} # if j is an headwater, impose StreamOrder = 1
        }
      } else if (streamOrderType=="Shreve"){
        # calculate Shreve stream order
        StreamOrder_AG <- numeric(vec_nNodesAG[thr])
        for (i in 1:vec_nNodesAG[thr]){
          j <- order(Nupstream_AG)[i] # index that explores reaches in a downstream direction
          tmp <- which(DownNode_AG==j) # set of reaches draining into j
          if (length(tmp)>0){
            StreamOrder_AG[j] <- sum(StreamOrder_AG[tmp])
          } else {StreamOrder_AG[j] <- 1} # if j is an headwater, impose StreamOrder = 1
        } 
      }
      vec_StreamOrder[thr] <- max(StreamOrder_AG)
    }
  }
  
  Thresholds <- vector("list",0)
  Thresholds[["thrValues"]] <- thrValues
  Thresholds[["nNodesRN"]] <- vec_nNodesRN
  Thresholds[["nNodesAG"]] <- vec_nNodesAG
  Thresholds[["drainageDensity"]] <- vec_DrainageDensity
  Thresholds[["streamOrder"]] <- vec_StreamOrder
  
  invisible(Thresholds)
}
