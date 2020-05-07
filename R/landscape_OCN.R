
landscape_OCN <- function(OCN,
                          slope0=1,
                          zMin=0,
                          optimizeDZ=FALSE,
                          optimMethod="BFGS",
                          optimControl=list(maxit=100*length(OCN$FD$outlet), trace=1),
                          displayUpdates=0) {
  
  if (!(displayUpdates %in% c(0,1,2))) {stop("Invalid displayUpdates")}
  
  if (displayUpdates>0){message("Calculating lengths and slopes... \r", appendLF = FALSE)}
  
  AvailableNodes <- setdiff(1:OCN$FD$nNodes,OCN$FD$outlet)  
  # calculate elevation gain through each pixel
  Slope <- slope0*(OCN$FD$A/(OCN$FD$nNodes*OCN$cellsize^2))^(OCN$expEnergy-1)
  Length <- rep(0,OCN$FD$nNodes)
  kount <- 0
  for (i in AvailableNodes){
    Length[i] <- sqrt((abs(OCN$FD$X[OCN$FD$downNode[i]]-OCN$FD$X[i]) %% ((OCN$dimX-1)*OCN$cellsize-2*min(OCN$FD$X)))^2 + 
                        (abs(OCN$FD$Y[OCN$FD$downNode[i]]-OCN$FD$Y[i]) %% ((OCN$dimY-1)*OCN$cellsize-2*min(OCN$FD$Y)))^2)
    kount <- kount + 1
    if (displayUpdates==2){message(sprintf("Calculating lengths and slopes... %.1f%%\r",kount/length(AvailableNodes)*100), appendLF = FALSE)}
  }
  DeltaZ <- Slope*Length
  
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
  if (displayUpdates>0){message("Calculating lengths and slopes...   100%\n", appendLF = FALSE)}
  
  # find elevation pattern with respect to main outlet
  if (displayUpdates>0){message("Determining elevation... \r", appendLF = FALSE)}
  kount <- 0
  Z <- numeric(OCN$FD$nNodes)
  FD_to_CM <- numeric(OCN$FD$nNodes)
  CM_to_FD <- vector("list",OCN$nOutlet)
  for (outlet in 1:length(OCN$FD$outlet)){
    next_nodes <- OCN$FD$outlet[outlet]
    FD_to_CM[OCN$FD$outlet[outlet]] <- outlet
    CM_to_FD[[outlet]] <- OCN$FD$outlet[outlet]
    while (length(next_nodes)>0) {
      current_nodes <- next_nodes
      kount <- kount + length(current_nodes)
      next_nodes <- integer(0) # empty next_nodes
      for (i in 1:length(current_nodes)){
        node <- current_nodes[i]
        neighbours <- which(OCN$FD$downNode==node)
        Z[neighbours] <- Z[node] + DeltaZ[neighbours]
        FD_to_CM[neighbours] <- outlet
        CM_to_FD[[outlet]] <- c(CM_to_FD[[outlet]],neighbours)
        next_nodes <- c(next_nodes,neighbours)
      }
      if (displayUpdates==2){message(sprintf("Determining elevation... %.1f%%\r",kount/OCN$FD$nNodes*100), appendLF = FALSE)}
    }
  }
  
  # determine catchment area
  A <- numeric(OCN$nOutlet)
  for (i in 1:OCN$nOutlet){
    A[i] <- sum(FD_to_CM==i)*OCN$cellsize^2
  }
  sortA <- sort(A,decreasing=TRUE,index.return=TRUE)
  
  # adjacency matrix at catchment level
  if (OCN$nOutlet>1){
    #W_CM <- sparseMatrix(i=1,j=1,x=0,dims=c(OCN$nOutlet,OCN$nOutlet))
    W_CM <- spam(0,OCN$nOutlet,OCN$nOutlet)
    for (i in 1:OCN$nOutlet){
      for (k in 1:length(CM_to_FD[[i]])){
        ind <- CM_to_FD[[i]][k]
        set <- NeighbouringNodes[[ind]]
        NeighSubcatch <- FD_to_CM[set]
        NeighSubcatch <- NeighSubcatch[!is.nan(NeighSubcatch)]
        Border <- which(NeighSubcatch!=i)
        if (length(Border)>0) {W_CM[i,unique(NeighSubcatch[Border])] <- 1}
      }
    }
  }else {W_CM <- 0}
  if (displayUpdates>0){message("Determining elevation...   100%\n", appendLF = FALSE)}
  
  
  # find altitude of secondary outlets with respect to altitude of the main outlet
  if (optimizeDZ==TRUE){
    if (displayUpdates==1){message("Optimizing outlet elevations... \r", appendLF = FALSE)}
    if (length(OCN$FD$outlet)>1){ 
      if (optimControl$trace>0) {message("Optimizing outlet elevations...\n", appendLF = FALSE)}
      CatchmentMat <- matrix(data=FD_to_CM,nrow=OCN$dimY,ncol=OCN$dimX)
      # find border pixels between catchments
      # BorderMat <- sparseMatrix(i=1,j=1,x=0,dims=c(OCN$FD$nNodes,OCN$FD$nNodes))
      BorderMat <- spam(0,OCN$FD$nNodes,OCN$FD$nNodes)
      ind <- matrix(0,1000*OCN$FD$nNodes,2)
      k <- 1
      for (i in 1:OCN$FD$nNodes){
        NeighCatch <- FD_to_CM[NeighbouringNodes[[i]]]
        isBorder <- (NeighCatch!=FD_to_CM[i])
        len <- length(NeighCatch[isBorder])
        if (len>0){
          # BorderMat[i,NeighbouringNodes[[i]][isBorder]] <- 1
          if ((k+len-1) <= dim(ind)[1]){
            ind[k:(k+len-1),] <- matrix(c(rep(i,len),NeighbouringNodes[[i]][isBorder]),nrow=len,ncol=2)
          } else {ind <- rbind(ind,matrix(c(rep(i,len),NeighbouringNodes[[i]][isBorder]),nrow=len,ncol=2))}
          
          k <- k + len 
          
        }
      }
      ind <- ind[-which(ind[,1]==0),]
      BorderMat[ind] <- 1
      
      # function for minimization of delta Z at the catchment borders
      OptimizeDeltaZ <- function(x) {
        UpdateZ <- rep(0,OCN$FD$nNodes)
        # the elevation of the biggest catchment is not changed
        for (i in 1:(length(OCN$FD$outlet)-1)){
          UpdateZ <- UpdateZ + (FD_to_CM==(sortA$ix[i+1]))*x[i]}
        Znew <- Z+UpdateZ
        # Znew <- Znew %*% t(1+numeric(length(Z)))
        mat <- BorderMat
        mat@entries <- mat@entries*rep(Znew, diff(mat@rowpointers))
        sum(abs(mat  - t(mat) )) # functional to be minimized
      }
      
      # use Nelder-Mead solver for minimization of OptimizeDeltaZ
      OptList <- optim(rep(0,length(OCN$FD$outlet)-1),OptimizeDeltaZ,method=optimMethod,
                       control=optimControl)
      Z_lifts <- OptList$par
      
      # apply lifting to catchments
      UpdateZ <- rep(0,OCN$FD$nNodes)
      for (i in 1:(length(OCN$FD$outlet)-1)){
        UpdateZ <- UpdateZ + (FD_to_CM==(sortA$ix[i+1]))*Z_lifts[i]}
      Z <- Z+UpdateZ
      
      if (min(Z_lifts)<0) {
        Z <- Z - min(Z_lifts)
        #print(sprintf("Outlet of main catchment has been lifted by %.2f elevation units",- min(Z_lifts)))
      }
    }
    if (displayUpdates>0){message("Optimizing outlet elevations... 100%\n", appendLF = FALSE)}
  }
  Z <- Z + zMin
  
  
  # exact drawing for reflecting boundaries networks
  X_draw <- OCN$FD$X # new vector of X coordinates
  Y_draw <- OCN$FD$Y # new vector of Y coordinates
  
  if(OCN$periodicBoundaries==TRUE){
    if (displayUpdates>0){message("Calculating real X, Y coordinates... \r", appendLF = FALSE)}
    kount <- 0
    CurrentPath <- OCN$FD$outlet # start reattributing coordinates from outlet(s)
    while (length(CurrentPath)>0){ # iterate until all pixels have been explored
      ContinuePath <- vector(mode="numeric", length=0) # create empty set of pixels upstream of CurrentPath
      for (k in 1:length(CurrentPath)){
        
        UpOneLevel <- which(OCN$FD$downNode==CurrentPath[k]) # find pixels upstream of CurrentPath
        if (length(UpOneLevel)>0){ # if CurrentPath[k] is not a headwater, continue
          for (i in 1:length(UpOneLevel)){
            # reattribute X coordinate of UpOneLevel[i]
            if (X_draw[UpOneLevel[i]]-X_draw[CurrentPath[k]] > OCN$cellsize){
              X_draw[UpOneLevel[i]] <- X_draw[UpOneLevel[i]] - 
                OCN$cellsize*OCN$dimX*(1 + floor((X_draw[UpOneLevel[i]] - X_draw[CurrentPath[k]] - 2*OCN$cellsize)/(OCN$dimX*OCN$cellsize))) 
              # the factor floor(...) is added to adjust for pixels that are flipped several times
            } else if (X_draw[UpOneLevel[i]]-X_draw[CurrentPath[k]] < -OCN$cellsize) {
              X_draw[UpOneLevel[i]] <- X_draw[UpOneLevel[i]] + 
                OCN$cellsize*OCN$dimX*(1+floor((X_draw[CurrentPath[k]] - X_draw[UpOneLevel[i]] - 2*OCN$cellsize)/(OCN$dimX*OCN$cellsize)))
            }
            # reattribute Y coordinate of UpOneLevel[i]
            if (Y_draw[UpOneLevel[i]]-Y_draw[CurrentPath[k]] > OCN$cellsize){
              Y_draw[UpOneLevel[i]] <- Y_draw[UpOneLevel[i]] - 
                OCN$cellsize*OCN$dimY*(1+floor((Y_draw[UpOneLevel[i]] - Y_draw[CurrentPath[k]] - 2*OCN$cellsize)/(OCN$dimY*OCN$cellsize)))
            } else if (Y_draw[UpOneLevel[i]]-Y_draw[CurrentPath[k]] < -OCN$cellsize) {
              Y_draw[UpOneLevel[i]] <- Y_draw[UpOneLevel[i]] + 
                OCN$cellsize*OCN$dimY*(1+floor((Y_draw[CurrentPath[k]] - Y_draw[UpOneLevel[i]] - 2*OCN$cellsize)/(OCN$dimY*OCN$cellsize)))
            }
          }
        }
        ContinuePath <- c(ContinuePath,UpOneLevel) # add UpOneLevel to set of pixels that will be explored on the next iteration
      }
      CurrentPath <- ContinuePath # move to next iteration
      kount <- kount + length(CurrentPath)
      if (displayUpdates==2){message(sprintf("Calculating real X, Y coordinates... %.1f%%\r",kount/OCN$FD$nNodes*100), appendLF = FALSE)}
    }
  }
  if (displayUpdates>0){message("Calculating real X, Y coordinates...   100%\n", appendLF = FALSE)}
  
  if (displayUpdates>0){message("Calculating catchment contour(s)... \r", appendLF = FALSE)}
  # determine contour of catchments (with original coordinates)
  X_contour <- vector("list", length = OCN$nOutlet)
  Y_contour <- vector("list", length = OCN$nOutlet)
  kount <- 0
  for (j in 1:OCN$nOutlet){
    subset_X <- OCN$FD$X[FD_to_CM==j]
    subset_Y <- OCN$FD$Y[FD_to_CM==j]
    X_mesh <- seq(min(subset_X)-OCN$cellsize,max(subset_X)+OCN$cellsize,OCN$cellsize)
    Y_mesh <- seq(min(subset_Y)-OCN$cellsize,max(subset_Y)+OCN$cellsize,OCN$cellsize)
    mesh <- matrix(data=0,nrow=length(Y_mesh),ncol=length(X_mesh))
    for (i in 1:OCN$FD$nNodes){
      ind_x <- which(X_mesh==subset_X[i])
      ind_y <- which(Y_mesh==subset_Y[i])
      mesh[ind_y,ind_x] <- 1
    }
    count <- contourLines(X_mesh,Y_mesh,t(mesh),levels=1)
    X_contour[[j]] <- vector("list", length = length(count))
    Y_contour[[j]] <- vector("list", length = length(count))
    for (k in 1:length(count)){
      X_contour[[j]][[k]] <- count[[k]]$x
      Y_contour[[j]][[k]] <- count[[k]]$y
    } 
    if (displayUpdates==2){message(sprintf("Calculating catchment contour(s)... %.1f%%\r",j/OCN$nOutlet*50), appendLF = FALSE)}
  }
  
  # determine contour of catchments (for real-shaped OCN)
  X_contour_draw <- vector("list", length = OCN$nOutlet)
  Y_contour_draw <- vector("list", length = OCN$nOutlet)
  for (j in 1:OCN$nOutlet){
    subset_X <- X_draw[FD_to_CM==j]
    subset_Y <- Y_draw[FD_to_CM==j]
    X_mesh <- seq(min(subset_X)-OCN$cellsize,max(subset_X)+OCN$cellsize,OCN$cellsize)
    Y_mesh <- seq(min(subset_Y)-OCN$cellsize,max(subset_Y)+OCN$cellsize,OCN$cellsize)
    mesh <- matrix(data=0,nrow=length(Y_mesh),ncol=length(X_mesh))
    for (i in 1:OCN$FD$nNodes){
      ind_x <- which(X_mesh==subset_X[i])
      ind_y <- which(Y_mesh==subset_Y[i])
      mesh[ind_y,ind_x] <- 1
    }
    count <- contourLines(X_mesh,Y_mesh,t(mesh),levels=1)
    X_contour_draw[[j]] <- vector("list", length = length(count))
    Y_contour_draw[[j]] <- vector("list", length = length(count))
    for (k in 1:length(count)){
      X_contour_draw[[j]][[k]] <- count[[k]]$x
      Y_contour_draw[[j]][[k]] <- count[[k]]$y
    }  
    if (displayUpdates==2){message(sprintf("Calculating catchment contour(s)... %.1f%%\r",50+j/OCN$nOutlet*50), appendLF = FALSE)}
  }
  if (displayUpdates==1){message("Calculating catchment contour(s)... 100%\n", appendLF = FALSE)}
  
  # add results to OCN list
  OCN$FD[["slope"]] <- Slope
  OCN$FD[["leng"]] <- Length
  OCN$FD[["toCM"]] <- FD_to_CM
  OCN$FD[["XDraw"]] <- X_draw
  OCN$FD[["YDraw"]] <- Y_draw
  OCN$FD[["Z"]] <- Z
  OCN$CM[["A"]] <- A
  OCN$CM[["W"]] <- W_CM
  OCN$CM[["XContour"]] <- X_contour
  OCN$CM[["YContour"]] <- Y_contour
  OCN$CM[["XContourDraw"]] <- X_contour_draw
  OCN$CM[["YContourDraw"]] <- Y_contour_draw
  
  
  OCN$slope0 <- slope0
  OCN$zMin <- zMin
  
  if (optimizeDZ==TRUE) {OCN$optList <- OptList}
  
  invisible(OCN)
}