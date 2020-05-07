
aggregate_OCN <- function(OCN,
                          thrA=0.002*OCN$dimX*OCN$dimY*OCN$cellsize^2,
                          streamOrderType="Strahler",
                          maxReachLength=Inf){
  
  if (!("slope" %in% names(OCN$FD))){
    stop('Missing fields in OCN. You should run landscape_OCN prior to aggregate_OCN.')
  }
  
  if (maxReachLength < OCN$cellsize*sqrt(2)){
    stop("maxReachLength cannot be smaller than OCN$cellsize*sqrt(2).")
  }
  #t1 <- Sys.time()
  
  ###############################
  ## BUILD NETWORK AT RN LEVEL ##
  ###############################
  
  #print('Crop data at FD level to RN level...',quote=FALSE); 
  RN_mask <- as.vector(OCN$FD$A >= thrA)# RN_mask allows to sample RN-level values from matrices/vectors at FD level   
  RN_to_FD <- which(RN_mask) # RN_to_FD[i] is the pixel ID at the FD level of the pixel whose ID at the RN level is i
  FD_to_RN <- RN_mask*cumsum(as.numeric(RN_mask)) # FD_to_RN[i] is the pixel ID at the RN level of the pixel whose ID at the FD level is i
  # if pixel i at FD level doesn't belong to RN, then FD_to_RN[i]=0
  
  Nnodes_RN <- length(RN_to_FD)
  
  W_RN <- OCN$FD$W[RN_mask,,drop=FALSE]
  W_RN <- W_RN[,RN_mask,drop=FALSE]
  
  Outlet_RN <- FD_to_RN[OCN$FD$outlet]
  Outlet_RN <- Outlet_RN[Outlet_RN!=0] # remove outlets if the corresponding catchment size is lower than thrAeshold
  DownNode_RN <- numeric(Nnodes_RN)
  # for (i in 1:Nnodes_RN){
  #   if (!(i %in% Outlet_RN)){
  #     DownNode_RN[i] <- which(W_RN[i,]==1)
  #   }}
  tmp <- W_RN@rowpointers
  NotOutlet <- which((tmp[-1] - tmp[-length(tmp)])==1)
  DownNode_RN[NotOutlet] <- W_RN@colindices
  
  A_RN <- OCN$FD$A[RN_mask]
  X_RN <- OCN$FD$X[RN_mask]
  Y_RN <- OCN$FD$Y[RN_mask]
  Z_RN <- OCN$FD$Z[RN_mask]
  Length_RN <- OCN$FD$leng[RN_mask]
  
  # Drainage density
  DrainageDensity_RN <- sum(Length_RN)/(OCN$dimX*OCN$dimY*OCN$cellsize^2)
  
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
  
  # Calculate slope for each pixel of the river network 
  Slope_RN <- OCN$FD$slope[RN_mask]
  #print(sprintf('Elapsed time %.2f s',difftime(Sys.time(),t1,units='secs')),quote=FALSE)
  #t1 <- Sys.time()
  
  # Upstream_RN : list containing IDs of all nodes upstream of each node (plus node itself)
  Upstream_RN <- vector("list",Nnodes_RN)
  Nupstream_RN <- numeric(Nnodes_RN)
  for (i in 1:Nnodes_RN){
    UpOneLevel <- which(DownNode_RN==i) # find reaches at one level upstream
    Upstream_RN[[i]] <- UpOneLevel      # add them to the list
    while (length(UpOneLevel)!=0) { # continue until there are no more reaches upstream
      ContinuePath <- UpOneLevel # jump 1 level above
      UpOneLevel <- which(DownNode_RN %in% ContinuePath) # find reaches at one level upstream
      Upstream_RN[[i]] <- c(Upstream_RN[[i]],UpOneLevel) # add them to the list
    }
    Upstream_RN[[i]] <- c(Upstream_RN[[i]],i)
    Nupstream_RN[i] <- length(Upstream_RN[[i]])
  }
  # RN_to_CM[i] indicates outlet to which reach i drains
  RN_to_CM <- numeric(Nnodes_RN)
  for (i in 1:OCN$nOutlet){
    RN_to_CM[Upstream_RN[[Outlet_RN[i]]]] <- i
  }
  
  
  ###############################
  ## BUILD NETWORK AT AG LEVEL ##
  ###############################
  
  # Vector that attributes reach ID to all river network pixels
  #print('Define nodes of aggregated network...',quote=FALSE); 
  Nnodes_AG <- sum(IsNodeAG)
  Length_AG <- numeric(Nnodes_AG)
  RN_to_AG <- numeric(Nnodes_RN)
  reachID <- 1
  X_AG <- NaN*numeric(Nnodes_AG)
  Y_AG <- NaN*numeric(Nnodes_AG)
  Z_AG <- NaN*numeric(Nnodes_AG)
  A_AG <- NaN*numeric(Nnodes_AG)
  while (length(whichNodeAG) != 0){ # explore all AG Nodes
    i <- whichNodeAG[1] # select the first
    RN_to_AG[i] <- reachID 
    j <- DownNode_RN[i] 
    X_AG[reachID] <- X_RN[i]
    Y_AG[reachID] <- Y_RN[i]
    Z_AG[reachID] <- Z_RN[i]
    A_AG[reachID] <- A_RN[i]
    Length_AG[reachID] <- Length_RN[i]
    tmp_length <- Length_RN[i]
    tmp <- NULL
    j0 <- j
    while (!IsNodeAG[j] && j!=0) {
      tmp <- c(tmp, j)
      tmp_length <-  tmp_length + Length_RN[j]
      j_old <- j
      j <- DownNode_RN[j]} 
    
    if (tmp_length > maxReachLength){
      n_splits <- ceiling(tmp_length/maxReachLength)
      new_maxLength <- tmp_length/n_splits
      j <- j0
      while (!IsNodeAG[j] && j!=0 && Length_AG[reachID] <= new_maxLength) {
        RN_to_AG[j] <- reachID 
        Length_AG[reachID] <-  Length_AG[reachID] + Length_RN[j]
        j_old <- j
        j <- DownNode_RN[j]}
      if (Length_AG[reachID] > new_maxLength){
        j <- j_old
        Length_AG[reachID] <-  Length_AG[reachID] - Length_RN[j]
        ChannelHeads[j] <- 1
        whichNodeAG <- c(whichNodeAG,j)}

    } else {
      RN_to_AG[tmp] <- reachID
      Length_AG[reachID] <- tmp_length
    }
    
    reachID <- reachID + 1
    whichNodeAG <- whichNodeAG[-1]
  }
  Nnodes_AG <- length(X_AG)
  #Nnodes_AG <- length(X_AG) + sum(OutletNotChannelHead) # recalculate number of nodes to account for new channel heads
  # if (sum(OutletNotChannelHead)>0){
  #   Length_AG[reachID:Nnodes_AG] <- 0
  # }
  # if thrA=0, do not perform aggregation. Every pixel is a node 
  # (note that with thrA=1, reaches with more than one pixel can exist)
  if (thrA==0){
    Nnodes_AG <- Nnodes_RN
    RN_to_AG <- 1:Nnodes_RN
  }
  
  # FD_to_SC: vector of length OCN$FD$nNodes containing subcatchmentID for every pixel of the catchment
  # AG_to_FD: list containing FD indices of pixels belonging to a given reach 
  # SC_to_FD: list containing FD indices of pixels belonging to a given subcatchment 
  FD_to_SC <- NaN*numeric(OCN$FD$nNodes)
  
  # initialize FD_to_SC by attributing SC values to pixels belonging to AG level
  FD_to_SC[RN_mask] <- RN_to_AG 
  
  # attribute new SC values to pixels corresponding to outlets of catchments without reaches (because the drained area of the catchment is < thrA)
  Nnodes_SC <- Nnodes_AG + sum(OCN$FD$A[OCN$FD$outlet]<thrA)
  FD_to_SC[OCN$FD$outlet[OCN$FD$A[OCN$FD$outlet] < thrA]] <- (Nnodes_AG+1):Nnodes_SC 
  IndexHeadpixel <- which(OCN$FD$A==OCN$cellsize^2) # find FD pixels corresponding to headwaters
  
  
  AG_to_FD <- vector("list", Nnodes_AG)
  AG_to_RN <- vector("list", Nnodes_AG)
  for(i in 1:Nnodes_AG) { # attribute river network pixels to fields of the AG_to_FD list 
    AG_to_FD[[i]] <- RN_to_FD[which(RN_to_AG==i)]
    AG_to_RN[[i]] <- which(RN_to_AG==i) 
  }
  SC_to_FD <- AG_to_FD[1:Nnodes_AG] # initialize SC_to_FD by attributing the pixels that belong to reaches
  
  # add pixels corresponding to outlets of catchments without reaches
  if (Nnodes_SC > Nnodes_AG){
    for (i in (Nnodes_AG+1):Nnodes_SC){
      SC_to_FD[[i]] <- OCN$FD$outlet[OCN$FD$A[OCN$FD$outlet]<thrA][i-Nnodes_AG]
    }}
  
  for (i in 1:length(IndexHeadpixel)){ # i: index that spans all headwater pixels
    p <- IndexHeadpixel[i] # p: ID of headwater pixel
    pNew <- p; # pNew: pixel downstream of p 
    k <- NaN; # k: SC value of pixel pNew
    sub_p <- integer(0) # sub_p is the subset of pixels downstream of pixel p
    while (is.nan(k)){ # continue downstream movement until a pixel to which the SC has already been attributed is found
      k <- FD_to_SC[pNew]
      if (is.nan(k)){
        sub_p <- c(sub_p,pNew)
        pNew <- OCN$FD$downNode[pNew]
      }}
    FD_to_SC[sub_p] <- k
    SC_to_FD[[k]] <- c(SC_to_FD[[k]],sub_p)
  }
  
  ######################################
  ## CALCULATE PROPERTIES AT AG LEVEl ##
  ######################################
  
  #print('W matrix at AG level...',quote=FALSE); 
  # Adjacency matrix at reach level
  DownNode_AG <- numeric(Nnodes_AG)
  # W_AG <- sparseMatrix(i=1,j=1,x=0,dims=c(Nnodes_AG,Nnodes_AG))
  W_AG <- spam(0,Nnodes_AG,Nnodes_AG)
  ind <- matrix(0,Nnodes_AG,2)
  reachID <- sum(ChannelHeads) + 1
  for (i in 1:Nnodes_RN){ 
    if (DownNode_RN[i] != 0 && RN_to_AG[DownNode_RN[i]] != RN_to_AG[i]) {
      DownNode_AG[RN_to_AG[i]] <- RN_to_AG[DownNode_RN[i]]
      #W_AG[RN_to_AG[i],DownNode_AG[RN_to_AG[i]]] <- 1
      ind[RN_to_AG[i],] <- c(RN_to_AG[i],DownNode_AG[RN_to_AG[i]])
    }
    # contributing area of nodes at AG level
    # if (ChannelHeads[i]){
    #   A_AG[RN_to_AG[i]] <- A_RN[i]
    # } 
  }
  ind <- ind[-which(ind[,1]==0),]
  W_AG[ind] <- 1
  Outlet_AG <- RN_to_AG[Outlet_RN]
  
  # Upstream_AG : list containing IDs of all reaches upstream of each reach (plus reach itself)
  Upstream_AG <- vector("list",Nnodes_AG)
  Nupstream_AG <- numeric(Nnodes_AG)
  for (i in 1:Nnodes_AG){
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
  # AG_to_CM[i] indicates outlet to which reach i drains
  AG_to_CM <- numeric(Nnodes_AG)
  for (i in 1:OCN$nOutlet){
    AG_to_CM[Upstream_AG[[Outlet_AG[i]]]] <- i
  }
  #print(sprintf('Elapsed time %.2f s',difftime(Sys.time(),t1,units='secs')),quote=FALSE)
  #t1 <- Sys.time()
  
  #print('Stream order at AG level...',quote=FALSE)
  if (streamOrderType=="Strahler"){
    # calculate Strahler stream order
    StreamOrder_AG <- numeric(Nnodes_AG)
    for (i in 1:Nnodes_AG){
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
    StreamOrder_AG <- numeric(Nnodes_AG)
    for (i in 1:Nnodes_AG){
      j <- order(Nupstream_AG)[i] # index that explores reaches in a downstream direction
      tmp <- which(DownNode_AG==j) # set of reaches draining into j
      if (length(tmp)>0){
        StreamOrder_AG[j] <- sum(StreamOrder_AG[tmp])
      } else {StreamOrder_AG[j] <- 1} # if j is an headwater, impose StreamOrder = 1
    } 
  }
  #print(sprintf('Elapsed time %.2f s',difftime(Sys.time(),t1,units='secs')),quote=FALSE); 
  #t1 <- Sys.time()
  
  #print('Length and slope at AG level...',quote=FALSE) 
  # Calculate length and slopes of reaches
  #Length_AG <- rep(0,Nnodes_AG)
  Slope_AG <- numeric(Nnodes_AG)
  for (i in 1:Nnodes_AG){
    #Length_AG[i] <- sum(OCN$FD$leng[AG_to_FD[[i]]])
    Slope_AG[i] <- (Slope_RN[RN_to_AG==i] %*% Length_RN[RN_to_AG==i])/Length_AG[i] # scalar product between vector of slopes and lengths of nodes at RN level belonging to reach i 
  }
  
  
  ######################################
  ## CALCULATE PROPERTIES AT SC LEVEL ##
  ######################################
  
  #print(sprintf('Elapsed time %.2f s',difftime(Sys.time(),t1,units='secs')),quote=FALSE) 
  #t1 <- Sys.time()
  
  #print('Subcatchment properties...',quote=FALSE) 
  # calculate subcatchment properties: Local Elevation, Local Drained Area, Upstream Area
  Z_SC <- numeric(Nnodes_SC)
  Alocal_SC <- numeric(Nnodes_SC)
  for (i in 1:Nnodes_SC) {
    Z_SC[i] <- mean(OCN$FD$Z[SC_to_FD[[i]]])
    Alocal_SC[i] <- length(SC_to_FD[[i]])*OCN$cellsize^2
  }
  # drained area at AG level: note that the first Nnodes_AG elements of Alocal_SC correspond to subcatchments with reaches
  
  #  Areach_AG: includes the areas drained by the reaches
  Areach_AG <- numeric(Nnodes_AG)
  for (i in 1:Nnodes_AG) {
    Areach_AG[i] <- sum(Alocal_SC[Upstream_AG[[i]]])  
  }

  # coordinates of AG nodes considered at the downstream end of the respective edge
  XReach <- numeric(Nnodes_AG)
  YReach <- numeric(Nnodes_AG)
  ZReach <- numeric(Nnodes_AG)
  for (i in 1:Nnodes_AG){
    tmp <- AG_to_RN[[i]]
    ind <- which(A_RN[tmp]==max(A_RN[tmp]))
    node <- tmp[ind]
    XReach[i] <- X_RN[node]
    YReach[i] <- Y_RN[node]
    ZReach[i] <- Z_RN[node]
  }
  XReach[Outlet_AG] <- NaN
  YReach[Outlet_AG] <- NaN
  ZReach[Outlet_AG] <- NaN
  
  # build neighbouring nodes at FD level
  # find list of possible neighbouring pixels
  movement <- matrix(c(0,-1,-1,-1,0,1,1,1,1,1,0,-1,-1,-1,0,1),nrow=2,byrow=TRUE)
  NeighbouringNodes <- vector("list", OCN$FD$nNodes)
  cont_node <- 0
  for (cc in 1:OCN$dimX) {
    for (rr in 1:OCN$dimY) {
      cont_node <- cont_node + 1
      neigh_r <- rep(rr,8)+movement[1,]
      neigh_c <- rep(cc,8)+movement[2,]
      if (OCN$periodicBoundaries == TRUE){
        neigh_r[neigh_r==0] <- OCN$dimY
        neigh_c[neigh_c==0] <- OCN$dimX
        neigh_r[neigh_r>OCN$dimY] <- 1
        neigh_c[neigh_c>OCN$dimX] <- 1
      }
      NotAboundary <- neigh_r>0 & neigh_r<=OCN$dimY & neigh_c>0 & neigh_c<=OCN$dimX # only effective when periodicBoundaries=FALSE
      NeighbouringNodes[[cont_node]] <- neigh_r[NotAboundary] + (neigh_c[NotAboundary]-1)*OCN$dimY
    }
  } 
  
  # Subcatchment adjacency matrix: find which subcatchments have borders in common
  #W_SC <- sparseMatrix(i=1,j=1,x=0,dims=c(Nnodes_SC,Nnodes_SC))
  W_SC <- spam(0,Nnodes_SC,Nnodes_SC)
  indices <- matrix(0,Nnodes_SC,2)
  for (i in 1:Nnodes_SC){
    for (k in 1:length(SC_to_FD[[i]])){
      ind <- SC_to_FD[[i]][k]
      if (length(ind)>0) {
        set <- NeighbouringNodes[[ind]]
        NeighSubcatch <- FD_to_SC[set]
        NeighSubcatch <- NeighSubcatch[!is.nan(NeighSubcatch)]
        Border <- which(NeighSubcatch!=i)
        if (length(Border)>0) {
          W_SC[i,unique(NeighSubcatch[Border])] <- 1
          }}
    }
  }
  
  # X,Y of subcatchment centroids
  X_SC <- numeric(Nnodes_SC)
  Y_SC <- numeric(Nnodes_SC)
  for (i in 1:Nnodes_SC){
    X_SC[i] <- mean(OCN$FD$X[SC_to_FD[[i]]])
    Y_SC[i] <- mean(OCN$FD$Y[SC_to_FD[[i]]])
  }
  
  
  ######################
  ## EXPORT VARIABLES ##
  ######################
  
  #FD level
  OCN$FD[["toRN"]] <- FD_to_RN
  OCN$FD[["toSC"]] <- FD_to_SC
  
  # RN level
  OCN$RN[["A"]] <- A_RN
  OCN$RN[["W"]] <- W_RN
  OCN$RN[["downNode"]] <- DownNode_RN
  OCN$RN[["drainageDensity"]] <- DrainageDensity_RN
  OCN$RN[["leng"]] <- Length_RN
  OCN$RN[["nNodes"]] <- Nnodes_RN
  OCN$RN[["nUpstream"]] <- Nupstream_RN
  OCN$RN[["outlet"]] <- Outlet_RN
  OCN$RN[["slope"]] <- Slope_RN
  OCN$RN[["toFD"]] <- RN_to_FD
  OCN$RN[["toAGReach"]] <- RN_to_AG
  OCN$RN[["toCM"]] <- RN_to_CM
  OCN$RN[["upstream"]] <- Upstream_RN
  OCN$RN[["X"]] <- X_RN
  OCN$RN[["Y"]] <- Y_RN
  OCN$RN[["Z"]] <- Z_RN
  
  # AG level
  OCN$AG[["A"]] <- A_AG
  OCN$AG[["AReach"]] <- Areach_AG
  OCN$AG[["W"]] <- W_AG
  OCN$AG[["downNode"]] <- DownNode_AG
  OCN$AG[["leng"]] <- Length_AG
  OCN$AG[["nNodes"]] <- Nnodes_AG
  OCN$AG[["nUpstream"]] <- Nupstream_AG
  OCN$AG[["outlet"]] <- Outlet_AG
  OCN$AG[["slope"]] <- Slope_AG
  OCN$AG[["streamOrder"]] <- StreamOrder_AG
  OCN$AG[["ReachToFD"]] <- AG_to_FD
  OCN$AG[["ReachToRN"]] <- AG_to_RN
  OCN$AG[["toCM"]] <- AG_to_CM
  OCN$AG[["upstream"]] <- Upstream_AG
  OCN$AG[["X"]] <- X_AG
  OCN$AG[["XReach"]] <- XReach
  OCN$AG[["Y"]] <- Y_AG
  OCN$AG[["YReach"]] <- YReach
  OCN$AG[["Z"]] <- Z_AG
  OCN$AG[["ZReach"]] <- ZReach
  
  # SC level
  OCN$SC[["ALocal"]] <- Alocal_SC
  OCN$SC[["W"]] <- W_SC
  OCN$SC[["nNodes"]] <- Nnodes_SC
  OCN$SC[["toFD"]] <- SC_to_FD
  OCN$SC[["X"]] <- X_SC
  OCN$SC[["Y"]] <- Y_SC  
  OCN$SC[["Z"]] <- Z_SC
  
  # other
  OCN$thrA <- thrA
  OCN$streamOrderType <- streamOrderType
  OCN$maxReachLength <- maxReachLength
  
  # re-define AG_to_RN, AG_to_FD, RN_to_AG considering AG nodes as pixels and not reaches
  AG_to_FDnode <- numeric(Nnodes_AG)
  AG_to_RNnode <- numeric(Nnodes_AG)
  for (i in 1:Nnodes_AG){
    tmpFD <- AG_to_FD[[i]]
    AG_to_FDnode[i] <- tmpFD[OCN$FD$A[tmpFD]==min(OCN$FD$A[tmpFD])]
    tmpRN <- AG_to_RN[[i]]
    AG_to_RNnode[i] <- tmpRN[OCN$RN$A[tmpRN]==min(OCN$RN$A[tmpRN])]
  }
  RN_to_AGnode <- numeric(Nnodes_RN)
  for (i in 1:Nnodes_AG){
    RN_to_AGnode[AG_to_RNnode[i]] <- i
  }
 
  OCN$RN[["toAG"]] <- RN_to_AGnode
  OCN$AG[["toFD"]] <- AG_to_FDnode
  OCN$AG[["toRN"]] <- AG_to_RNnode
  
  invisible(OCN)
}