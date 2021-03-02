continue_OCN <- function(OCN,nNewIter,
                         coolingRate=NULL,
                         initialNoCoolingPhase=0,
                         displayUpdates=1,
                         showIntermediatePlots=FALSE,
                         thrADraw=NULL,
                         easyDraw=NULL,
                         nUpdates=50){
  
  
  nNodes <- OCN$dimX*OCN$dimY
  dim <- sqrt(nNodes)
  A <- OCN$FD$A/OCN$cellsize^2
  X <- OCN$FD$X
  Y <- OCN$FD$Y
  dimX <- OCN$dimX
  dimY <- OCN$dimY
  downNode <- OCN$FD$downNode
  cellsize <- OCN$cellsize
  nOutlet <- OCN$nOutlet
  OutletPixel <- OCN$FD$outlet
  periodicBoundaries <- OCN$periodicBoundaries
  expEnergy <- OCN$expEnergy
  W <- OCN$FD$W
  Wt <- t(W)
  
  if (is.null(coolingRate)){
    coolingRate <- OCN$coolingRate[length(OCN$coolingRate)]
    resetCoolingSchedule <- FALSE
  } else {resetCoolingSchedule <- TRUE}
  
  if (dim^2 > 1000) {
    estTime <- (7.988e-1 - 1.266e-1*dim + 7.198e-3*dim^2 - 6.807e-5*dim^3 +  1.372e-6*dim^4) / (30*dim^2) * nNewIter
    unitTime <- "seconds"
    if (estTime >= 60 && estTime < 3600 ){
      estTime <- estTime/60
      unitTime <- "minutes"
    } else if (estTime >= 3600) {
      estTime <- estTime/3600
      unitTime <- "hours"}
    
    if (displayUpdates != 0){
      message("continue_OCN is running...\n", appendLF = FALSE)
      message(sprintf("Estimated duration: %.2f %s \n",estTime,unitTime), appendLF = FALSE)
      message("Note that the above estimate is only based on the choice of parameter nNewIter, and not on processor performance.\n", appendLF = FALSE)
      message("\n", appendLF = FALSE)
    }
  }
  
  t0 <- Sys.time()
  if (displayUpdates==2){message('Initializing...\n', appendLF = FALSE)}
  
  if (is.null(easyDraw)){
    if (dimX*dimY>4e4) {
      easyDraw=TRUE
    } else {easyDraw=FALSE}
  }
  
  if (is.null(thrADraw)){
    thrADraw <- 0.002*OCN$dimX*OCN$dimY*OCN$cellsize^2
  }
  
  if (is.null(OCN$FD$perm)){
    message("OCN$FD$perm is missing. Recalculating...\n", appendLF = FALSE)
    message("\n", appendLF = FALSE)
    pl <- initial_permutation(downNode)
    pl <- pl$perm
  } else {pl <- OCN$FD$perm}
  
  if(is.null(OCN$energyInit)){
    message("OCN$energyInit is missing. Recalculating...\n", appendLF = FALSE)
    message("\n", appendLF = FALSE)
    
    outletRow <- OutletPixel %% dimY
    outletRow[outletRow==0] <- dimY
    outletCol <- 1 + (OutletPixel-outletRow)/dimY
    
    outletPos <- numeric(nOutlet)
    outletSide <- character(nOutlet)
    for (i in 1:nOutlet){
      if (outletCol[i]==1){
        outletSide[i] <- "W"
        outletPos[i] <- outletRow[i]
      } else if (outletCol[i]==dimX){
        outletSide[i] <- "E"
        outletPos[i] <- outletRow[i]
      }
      if (outletRow[i]==1){
        outletSide[i] <- "S"
        outletPos[i] <- outletCol[i]
      } else if (outletRow[i]==dimY){
        outletSide[i] <- "N"
        outletPos[i] <- outletCol[i]
      }
    }
    
    ocn <- create_OCN(dimX,dimY,nIter=0,
                      outletPos=outletPos,outletSide=outletSide,nOutlet=OCN$nOutlet,typeInitialState=OCN$typeInitialState,displayUpdates=0)
    energyInit <- sum(ocn$FD$A^0.5)
  } else {energyInit <- OCN$energyInit}
  
  perimeter <- c(1:dimY, dimY*c(2:(dimX-1)), (dimX-1)*dimY+c(dimY:1), dimY*c((dimX-2):1)+1)
  perimeter <- rep(perimeter, 3)
  
  semiOutlet <- NULL
  if (dimX*dimY >= 1000){
    for (i in 1:nOutlet){
      tmp <- which(perimeter[(length(perimeter)/3+1):(length(perimeter)*2/3)] == OutletPixel[i])
      semiOutlet <- c(semiOutlet, perimeter[tmp - 1], perimeter[tmp + 1])  
    }
  }
  
  AvailableNodes <- setdiff(1:nNodes, union(OutletPixel,semiOutlet))
  AvailableNodesPlot <- setdiff(1:nNodes, OutletPixel)
  
  movement <- matrix(c(0,-1,-1,-1,0,1,1,1,1,1,0,-1,-1,-1,0,1),nrow=2,byrow=TRUE)
  NeighbouringNodes <- vector("list", nNodes)
  cont_node <- 0
  for (cc in 1:dimX) {
    for (rr in 1:dimY) {
      cont_node <- cont_node + 1
      neigh_r <- rr+movement[1,]
      neigh_c <- cc+movement[2,]
      if (periodicBoundaries==TRUE){
        neigh_r[neigh_r==0] <- dimY
        neigh_c[neigh_c==0] <- dimX
        neigh_r[neigh_r>dimY] <- 1
        neigh_c[neigh_c>dimX] <- 1
      }
      NotAboundary <- neigh_r>0 & neigh_r<=dimY & neigh_c>0 & neigh_c<=dimX # only effective when periodicBoundaries=FALSE
      NotAboundary_N4 <- neigh_r>0 & neigh_r<=dimY & neigh_c>0 & neigh_c<=dimX & (((1:8) %% 2)==1)
      NeighbouringNodes[[cont_node]] <- neigh_r[NotAboundary] + (neigh_c[NotAboundary]-1)*dimY
    }
  } 
  
  # plot initial state
  if (showIntermediatePlots==TRUE){
    AA <- A*cellsize^2
    if (nOutlet > 1){
      rnbw <- hcl.colors(nOutlet,palette="Dark 3")
      rnbw <- c(rnbw[(round(nOutlet/2)+1):nOutlet],rnbw[1:(round(nOutlet/2)+1)])
    } else {rnbw <- hcl.colors(3,palette="Dark 3")
    rnbw <- rnbw[3]}
    resort <- invPerm(pl)
    catch <- numeric(nNodes)
    for (o in 1:nOutlet){
      catch[(resort[OutletPixel[o]]-AA[OutletPixel[o]]/cellsize^2+1):resort[OutletPixel[o]]] <- o
    }
    old.par <- par(bty="n")
    on.exit(par(old.par))
    plot(c(min(X),max(X)),c(min(Y),max(Y)),main=sprintf('OCN %dx%d (initial state)',dimX,dimY),
         type="n",asp=1,axes=FALSE,xlab="",ylab="") # 
    points(X[OutletPixel],Y[OutletPixel],pch=15,col=rnbw[catch[resort[OutletPixel]]])
    if (easyDraw == FALSE){
      for (i in AvailableNodesPlot){
        if (AA[i]<=thrADraw & abs(X[i]-X[downNode[i]])<=cellsize & abs(Y[i]-Y[downNode[i]])<=cellsize ) {
          lines(c(X[i],X[downNode[i]]),c(Y[i],Y[downNode[i]]),lwd=0.5,col="#E0E0E0")}
      }
    }
    for (i in 1:nNodes){
      if (!(i %in% OutletPixel)){
        if (AA[i]>thrADraw & abs(X[i]-X[downNode[i]])<=cellsize & abs(Y[i]-Y[downNode[i]])<=cellsize ) {
          lines(c(X[i],X[downNode[i]]),c(Y[i],Y[downNode[i]]),lwd=0.5+4.5*(AA[i]/nNodes/cellsize^2)^0.5,col=rnbw[catch[resort[i]]])}
      }}
  }
  
  # initialize energy
  Energy <- numeric(nNewIter)
  Energy_0 <- sum(A^0.5)
  Energy[1] <-Energy_0
  savetime <- numeric(nNewIter)
  ExitFlag <- numeric(nNewIter)
  
  if (resetCoolingSchedule){
  Temperature <- c(energyInit + numeric(initialNoCoolingPhase*nNewIter),
                   energyInit*exp(-coolingRate*(1:(nNewIter-initialNoCoolingPhase*nNewIter))/nNodes))
  } else {
    Temperature <- c(energyInit + numeric(initialNoCoolingPhase*nNewIter),
                     energyInit*exp(-coolingRate*(OCN$nIter + (1:(nNewIter-initialNoCoolingPhase*nNewIter)))/nNodes))  
  }

  
  if (displayUpdates == 2){
    message(sprintf('Initialization completed. Elapsed time is %.2f s \n',difftime(Sys.time(),t0,units='secs')), appendLF = FALSE) 
    message('Search algorithm has started.\n', appendLF = FALSE)
    message('\n', appendLF = FALSE)}
  
  # simulated annealing algorithm
  if (nNewIter > 1){
    for (iter in 2:nNewIter){
      t2 <- Sys.time() 
      # pick random node (excluding the outlet)
      Energy_new <- 100*energyInit
      node <- sample(AvailableNodes,1)
      # change downstream connection from chosen node
      down_new <- sample(NeighbouringNodes[[node]][NeighbouringNodes[[node]]!=downNode[node]],1) # sample one node from list of neighbouring nodes, excluding the one that was previously connected to node
      
      ind1 <- (X[NeighbouringNodes[[node]]]==X[node] & Y[NeighbouringNodes[[node]]]==Y[down_new])
      Node1 <- NeighbouringNodes[[node]][ind1]
      
      ind2 <- (X[NeighbouringNodes[[node]]]==X[down_new] & Y[NeighbouringNodes[[node]]]==Y[node])
      Node2 <- NeighbouringNodes[[node]][ind2]
      
      if (isTRUE(Node1 != down_new) && isTRUE(Node2 != down_new) && # if  node -> down_new is diagonal connection 
          (downNode[Node1] == Node2 || downNode[Node2] == Node1)) { # if cross flow
        flag <- 3 # skip iteration because of cross-flow
      } else {
        
        pas <- allinoneF( as.integer(nOutlet), pl, Wt, downNode, node,
                          down_new, A[node], expEnergy)
        flag <- pas$flag
        if (flag==1){
          Anew <- pas$Anew
          Energy_new <- pas$energy
        }
      }
      
      # accept change if energy is lower or owing to the simulated annealing rule
      if (Energy_new <= Energy[iter-1] | runif(1)<exp(-(Energy_new-Energy[iter-1])/Temperature[iter-1])) {
        # update network
        flag <- 0
        pl <- pas$perm
        Wt <- pas$Wt_new
        A <- Anew[invPerm(pl)] # re-permute to the original indexing
        Energy[iter] <- Energy_new
        downNode[node] <- down_new
      } else {Energy[iter] <- Energy[iter-1]} # recject change and keep previous W
      
      # write update
      if (displayUpdates==2){
        if (iter %% round(nNewIter/nUpdates)==0){
          message(sprintf('%.1f%% completed - Elapsed time: %.2f s - %11s - Energy: %0.f \n',
                          iter/nNewIter*100,difftime(Sys.time(),t0,units='secs'),format(Sys.time(),"%b%d %H:%M"),Energy[iter]), appendLF = FALSE)
        }}
      # plot update
      
      if (iter %% round(nNewIter/nUpdates)==0){
        if (showIntermediatePlots==TRUE){ 
          AA <- A*cellsize^2 
          resort <- invPerm(pl)
          catch <- numeric(nNodes)
          for (o in 1:nOutlet){
            catch[(resort[OutletPixel[o]]-AA[OutletPixel[o]]/cellsize^2+1):resort[OutletPixel[o]]] <- o
          }
          plot(c(min(X),max(X)),c(min(Y),max(Y)),type="n",main=sprintf('OCN %dx%d (%.1f%% completed)',dimX,dimY,iter/nNewIter*100),
               asp=1,axes=FALSE,xlab=" ",ylab=" ") # 
          points(X[OutletPixel],Y[OutletPixel],pch=15,col=rnbw[catch[resort[OutletPixel]]])
          if (easyDraw == FALSE){
            for (i in AvailableNodesPlot){
              if (AA[i]<=(thrADraw)  & abs(X[i]-X[downNode[i]])<=cellsize & abs(Y[i]-Y[downNode[i]])<=cellsize ) {
                lines(c(X[i],X[downNode[i]]),c(Y[i],Y[downNode[i]]),lwd=0.5,col="#E0E0E0")}
            }
          }
          for (i in 1:nNodes){
            if (!(i %in% OutletPixel)){
              if (AA[i]>(thrADraw)  & abs(X[i]-X[downNode[i]])<=cellsize & abs(Y[i]-Y[downNode[i]])<=cellsize ) {
                lines(c(X[i],X[downNode[i]]),c(Y[i],Y[downNode[i]]),lwd=0.5+4.5*(AA[i]/nNodes/cellsize^2)^0.5,col=rnbw[catch[resort[i]]])}
            }}
        }}
      #}
      t3<-Sys.time()
      savetime[iter] <- t3-t2
      ExitFlag[iter] <- flag
      
      if (sum(A[OutletPixel])  != dimX*dimY){
        stop('Error: sum(A[OutletPixel]) is not equal to the total lattice area')
      }
    }
  }
  
  
  #W <- as.dgCMatrix.spam(t(Wt)) # ensure compatibility with other functions
  W <- t(Wt)
  
  FD <- list(A=A*cellsize^2,W=W,downNode=downNode,X=X,Y=Y,nNodes=nNodes,outlet=OutletPixel,perm=pl)
  OCN_new <- list(FD=FD,dimX=dimX,dimY=dimY,cellsize=cellsize,nOutlet=nOutlet,periodicBoundaries=periodicBoundaries,
                  expEnergy=expEnergy,coolingRate=c(OCN$coolingRate,coolingRate),typeInitialState=OCN$typeInitialState,nIter=OCN$nIter+nNewIter,
                  initialNoCoolingPhase=c(OCN$initialNoCoolingPhase,initialNoCoolingPhase),
                  energyInit=energyInit)
  
  if (!is.null(OCN$energy)) {OCN_new[["energy"]] <- c(OCN$energy,Energy)}
  if (!is.null(OCN$exitFlag)) {OCN_new[["exitFlag"]] <- c(OCN$exitFlag,ExitFlag)}
  if (!(is.null(OCN$N8))) {
    OCN_new[["N8"]] <- OCN$N8}
  if (!(is.null(OCN$N4))) {
    OCN_new[["N4"]] <- OCN$N4}
  
  invisible(OCN_new)
  
}