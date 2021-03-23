create_OCN <- function(dimX,dimY,
                       nOutlet=1,
                       outletSide="S", # vector of sides where outlets are located
                       outletPos=round(dimX/3),
                       periodicBoundaries=FALSE, # if TRUE, reflecting boundaries are applied. If FALSE, boundaries are non-reflecting.
                       typeInitialState=NULL,
                       flowDirStart=NULL,
                       expEnergy=0.5, # energy \propto Q*deltaH, Q propto A, deltaH propto A^-0.5
                       cellsize=1,
                       xllcorner=0.5*cellsize,
                       yllcorner=0.5*cellsize,
                       nIter=40*dimX*dimY,
                       nUpdates=50, # number of times an update is shown
                       initialNoCoolingPhase=0,
                       coolingRate=1,
                       showIntermediatePlots=FALSE,
                       thrADraw=0.002*dimX*dimY*cellsize^2,
                       easyDraw=NULL,
                       saveEnergy=FALSE,
                       saveExitFlag=FALSE,
                       saveN8=FALSE,
                       saveN4=FALSE,
                       displayUpdates=1){
  
  if (dimX<2 | dimY<2) stop("Dimensions too small.")
  if (dimX*dimY > 1000) {
    dim <- sqrt(dimX*dimY)
    estTime <- (7.988e-1 - 1.266e-1*dim + 7.198e-3*dim^2 - 6.807e-5*dim^3 +  1.372e-6*dim^4) / (30*dimX*dimY) * nIter
    unitTime <- "seconds"
    if (estTime >= 60 && estTime < 3600 ){
      estTime <- estTime/60
      unitTime <- "minutes"
    } else if (estTime >= 3600) {
      estTime <- estTime/3600
      unitTime <- "hours"}
    
    if (displayUpdates != 0){
      message("create_OCN is running...\n", appendLF = FALSE)
      message(sprintf("Estimated duration: %.2f %s \n",estTime,unitTime), appendLF = FALSE)
      message("Note that the above estimate is only based on the choice of parameters dimX, dimY and nIter, and not on processor performance.\n", appendLF = FALSE)
      message("\n", appendLF = FALSE)
    }
  }
  
  t0 <- Sys.time()
  if (displayUpdates==2){message('Initializing...\n', appendLF = FALSE)}
  
  ################################## 
  ## DEFINE INITIAL NETWORK STATE ##  
  ##################################  
  
  if ((nOutlet=="All")==TRUE && is.null(typeInitialState)==TRUE){
    typeInitialState <- "H"
  } else {
    if (is.null(typeInitialState)==TRUE){
      typeInitialState <- "I"
    }
  }
  
  # set outletSide, outletPos in accordance with the input given 
  if ((nOutlet=="All")==TRUE){
    outletSide <- c(rep("S",dimX),rep("E",dimY-2),rep("N",dimX),rep("W",dimY-2))
    outletPos <- c(1:dimX,2:(dimY-1),seq(dimX,1,-1),seq((dimY-1),2,-1))
    nOutlet <- 2*(dimX+dimY-2)
  } else if (nOutlet>1 && length(outletSide)==1 && length(outletPos)==1 && all(outletSide=="S") && all(outletPos==round(dimX/3))) {
    AllSides <- c(rep("S",dimX),rep("E",dimY-2),rep("N",dimX),rep("W",dimY-2))
    AllPos <- c(1:dimX,2:(dimY-1),1:dimX,2:(dimY-1))
    BorderPixels <- 1:(2*(dimX+dimY-2))
    outlet_ID <- sample(BorderPixels,nOutlet)
    outletSide <- AllSides[outlet_ID]
    outletPos <- AllPos[outlet_ID]
    # for (i in 1:nOutlet){
    #   outletSide[i] <- sample(AllSides,1)
    #   if (outletSide[i]=="N" || outletSide[i]=="S"){
    #     outletPos[i] <- sample(1:dimX,1) 
    #   } else {outletPos[i] <- sample(1:dimY,1)}
    # }
  }
  
  # check whether input is correct
  if(nOutlet != length(outletSide) || nOutlet != length(outletPos) ){
    stop('Length of outletSide and/or outletPos is inconsistent with nOutlet')}
  
  if (!all(outletSide %in% c("N","W","S","E"))){
    stop('Invalid outletSide')}
  
  for (i in 1:nOutlet){
    if (outletSide[i] %in% c("N","S")){
      if (!(outletPos[i] %in% 1:dimX)){
        stop('Invalid outletPos')} 
    } else if (!(outletPos[i] %in% 1:dimY)){
      stop('Invalid outletPos')}
  }
  
  if (!(typeInitialState %in% c("I","T","V","H"))){
    stop('Invalid typeInitialState')}
  
  if (is.null(easyDraw)){
    if (dimX*dimY>4e4) {
      easyDraw=TRUE
    } else {easyDraw=FALSE}
  }
  
  # define domain coordinates
  X <- rep(seq(xllcorner,xllcorner+(dimX-1)*cellsize,cellsize), each = dimY)  
  Y <- rep(seq(yllcorner,yllcorner+(dimY-1)*cellsize,cellsize),dimX)
  
  # find list of possible neighbouring pixels
  Nnodes <- length(X)
  movement <- matrix(c(0,-1,-1,-1,0,1,1,1,1,1,0,-1,-1,-1,0,1),nrow=2,byrow=TRUE)
  NeighbouringNodes <- vector("list", Nnodes)
  if (saveN8) {W_N8 <- spam(0, Nnodes,Nnodes)}  # null matrices
  if (saveN4) {W_N4 <- spam(0, Nnodes,Nnodes)}  
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
      if (saveN8){W_N8[cont_node,neigh_r[NotAboundary] + (neigh_c[NotAboundary]-1)*dimY]=1}
      if (saveN4){W_N4[cont_node,neigh_r[NotAboundary_N4] + (neigh_c[NotAboundary_N4]-1)*dimY]=1}
    }
  } 
  
  # find outlet position
  Outlet_row <- numeric(nOutlet)
  Outlet_col <- numeric(nOutlet)
  for (i in 1:nOutlet){
    if (outletSide[i]=="N"){
      Outlet_row[i] <- dimY
      Outlet_col[i] <- outletPos[i]
    } else if (outletSide[i]=="E") {
      Outlet_row[i] <- outletPos[i]
      Outlet_col[i] <- dimX
    } else if (outletSide[i]=="S") {
      Outlet_row[i] <- 1
      Outlet_col[i] <- outletPos[i]
    } else if (outletSide[i]=="W") {
      Outlet_row[i] <- outletPos[i]
      Outlet_col[i] <- 1
    }
  }
  
  # find flow direction matrix given initial network state
  if (is.null(flowDirStart)==TRUE){
    flowDirStart <- initialstate_OCN(dimX,dimY,nOutlet,outletSide,outletPos,typeInitialState)}
  
  # attribute flow direction =0 to outlet pixels
  for (i in 1:nOutlet){
    if (outletSide[i]=="N"){
      flowDirStart[dimY,outletPos[i]] <- 0 
      Outlet_row[i] <- dimY
      Outlet_col[i] <- outletPos[i]
    } else if (outletSide[i]=="E") {
      flowDirStart[outletPos[i],dimX] <- 0 
      Outlet_row[i] <- outletPos[i]
      Outlet_col[i] <- dimX
    } else if (outletSide[i]=="S") {
      flowDirStart[1,outletPos[i]] <- 0
      Outlet_row[i] <- 1
      Outlet_col[i] <- outletPos[i]
    } else if (outletSide[i]=="W") {
      flowDirStart[outletPos[i],1] <- 0
      Outlet_row[i] <- outletPos[i]
      Outlet_col[i] <- 1
    }
  }
  
  # Adjacency matrix for the initial state
  cont_node <- 0
  W <- spam(0,Nnodes,Nnodes)
  ind <- matrix(0,Nnodes,2)
  DownNode <- numeric(Nnodes)   
  for (cc in 1:dimX) {
    for (rr in 1:dimY) {
      cont_node <- cont_node + 1
      dir <- flowDirStart[rr,cc]
      if (dir>0) {
        DownPixel <- c(rr,cc)+movement[,dir]
        # if a custom flowDirStart is provided, allow flow to cross boundaries
        if (periodicBoundaries==TRUE){
          if (DownPixel[1]==0){DownPixel[1] <- dimY}
          if (DownPixel[1]>dimY){DownPixel[1] <- 1}
          if (DownPixel[2]==0){DownPixel[2] <- dimX}
          if (DownPixel[2]>dimX){DownPixel[2] <- 1}
        }
        ind[cont_node, ] <- c(cont_node,(DownPixel[2]-1)*dimY+DownPixel[1])
        #W[cont_node,(DownPixel[2]-1)*dimY+DownPixel[1]] <- 1
        DownNode[cont_node] <- (DownPixel[2]-1)*dimY+DownPixel[1]
      }
    }
  }
  ind <- ind[-which(ind[,1]==0),]
  W[ind] <- 1
  
  # patch to correct for initial 0
  newW <- new("spam")
  slot(newW, "entries", check = FALSE) <- W@entries[-1]
  slot(newW, "colindices", check = FALSE) <- W@colindices[-1]
  slot(newW, "rowpointers", check = FALSE) <- c(W@rowpointers[1],W@rowpointers[-1]-W@rowpointers[1])
  slot(newW, "dimension", check = FALSE) <- W@dimension 
  
  Wt <- t(newW)
  rm(W,newW)
  
  OutletPixel <- (Outlet_col-1)*dimY+Outlet_row # this doesn't keep the order
  
  perimeter <- c(1:dimY, dimY*c(2:(dimX-1)), (dimX-1)*dimY+c(dimY:1), dimY*c((dimX-2):1)+1)
  perimeter <- rep(perimeter, 3)
  
  semiOutlet <- NULL
  if (dimX*dimY >= 1000){
    for (i in 1:nOutlet){
      tmp <- which(perimeter[(length(perimeter)/3+1):(length(perimeter)*2/3)] == OutletPixel[i])
      semiOutlet <- c(semiOutlet, perimeter[tmp - 1], perimeter[tmp + 1])  
    }
  }
  
  AvailableNodes <- setdiff(1:Nnodes, union(OutletPixel,semiOutlet))
  AvailableNodesPlot <- setdiff(1:Nnodes, OutletPixel)
  
  pl <- initial_permutation(DownNode)  # calculate permutation vector
  
  pas <- permuteAddSolve(Wt, pl$perm, 1L, numeric(Nnodes), expEnergy)
  A <- pas[[1]]
  A <- A[invPerm(pl$perm)]
  
  # define identity matrix and vector of ones (needed to calculate A)
  #Imat <- diag.spam(1,Nnodes,Nnodes)
  #Ones <- numeric(1,Nnodes)
  # calculate vector of contributing area
  #A <- solve.spam(Imat-Wt,Ones)#*cellsize^2
  
  
  ## spam solution
  #Wt <- as.spam.dgCMatrix(Wt)
  
  
  ##########################
  ## OCN SEARCH ALGORITHM ##  
  ##########################  
  
  # initialize energy and temperatures
  Energy <- numeric(nIter)
  Energy_0 <- pas[[2]]
  Energy[1] <-Energy_0
  Temperature <- c(Energy[1]+numeric(initialNoCoolingPhase*nIter),Energy[1]*exp(-coolingRate*(1:(nIter-initialNoCoolingPhase*nIter))/(dimX*dimY)))
  savetime <- numeric(nIter)
  ExitFlag <- numeric(nIter)
  
  # plot initial state
  if (showIntermediatePlots==TRUE){
    AA <- A*cellsize^2
    if (nOutlet > 1){
      rnbw <- hcl.colors(nOutlet,palette="Dark 3")
      rnbw <- c(rnbw[(round(nOutlet/2)+1):nOutlet],rnbw[1:(round(nOutlet/2)+1)])
    } else {rnbw <- hcl.colors(3,palette="Dark 3")
    rnbw <- rnbw[3]}
    resort <- invPerm(pl$perm)
    catch <- numeric(Nnodes)
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
        if (AA[i]<=thrADraw & abs(X[i]-X[DownNode[i]])<=cellsize & abs(Y[i]-Y[DownNode[i]])<=cellsize ) {
          lines(c(X[i],X[DownNode[i]]),c(Y[i],Y[DownNode[i]]),lwd=0.5,col="#E0E0E0")}
      }
    }
    for (i in 1:Nnodes){
      if (!(i %in% OutletPixel)){
        if (AA[i]>thrADraw & abs(X[i]-X[DownNode[i]])<=cellsize & abs(Y[i]-Y[DownNode[i]])<=cellsize ) {
          lines(c(X[i],X[DownNode[i]]),c(Y[i],Y[DownNode[i]]),lwd=0.5+4.5*(AA[i]/Nnodes/cellsize^2)^0.5,col=rnbw[catch[resort[i]]])}
      }}
  }
  
  #Rprof("\\\\eawag/userdata/carrarlu/Desktop/Luca/OCN/R/spam_stuff/Rprof100.out")
  if (displayUpdates == 2){
    message(sprintf('Initialization completed. Elapsed time is %.2f s \n',difftime(Sys.time(),t0,units='secs')), appendLF = FALSE) 
    message('Search algorithm has started.\n', appendLF = FALSE)
    message('\n', appendLF = FALSE)}
  
  pl <- as.integer(pl$perm)
  flag <- 0
  
  # simulated annealing algorithm
  if (nIter > 1){
    for (iter in 2:nIter) {
      t2 <- Sys.time() 
      # pick random node (excluding the outlet)
      Energy_new <- 100*Energy_0
      node <- sample(AvailableNodes,1)
      # change downstream connection from chosen node
      down_new <- sample(NeighbouringNodes[[node]][NeighbouringNodes[[node]]!=DownNode[node]],1) # sample one node from list of neighbouring nodes, excluding the one that was previously connected to node
      
      ind1 <- (X[NeighbouringNodes[[node]]]==X[node] & Y[NeighbouringNodes[[node]]]==Y[down_new])
      Node1 <- NeighbouringNodes[[node]][ind1]
      
      ind2 <- (X[NeighbouringNodes[[node]]]==X[down_new] & Y[NeighbouringNodes[[node]]]==Y[node])
      Node2 <- NeighbouringNodes[[node]][ind2]
      
      if (isTRUE(Node1 != down_new) && isTRUE(Node2 != down_new) && # if  node -> down_new is diagonal connection 
          (DownNode[Node1] == Node2 || DownNode[Node2] == Node1)) { # if cross flow
        flag <- 3 # skip iteration because of cross-flow
      } else {
        
        pas <- allinoneF( as.integer(nOutlet), pl, Wt, DownNode, node,
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
        DownNode[node] <- down_new
      } else {Energy[iter] <- Energy[iter-1]} # recject change and keep previous W
      
      # write update
      if (displayUpdates==2){
        if (iter %% round(nIter/nUpdates)==0){
          message(sprintf('%.1f%% completed - Elapsed time: %.2f s - %11s - Energy: %0.f \n',
                          iter/nIter*100,difftime(Sys.time(),t0,units='secs'),format(Sys.time(),"%b%d %H:%M"),Energy[iter]), appendLF = FALSE)
        }}
      # plot update
      
      if (iter %% round(nIter/nUpdates)==0){
        if (showIntermediatePlots==TRUE){ 
          AA <- A*cellsize^2 
          resort <- invPerm(pl)
          catch <- numeric(Nnodes)
          for (o in 1:nOutlet){
            catch[(resort[OutletPixel[o]]-AA[OutletPixel[o]]/cellsize^2+1):resort[OutletPixel[o]]] <- o
          }
          plot(c(min(X),max(X)),c(min(Y),max(Y)),type="n",main=sprintf('OCN %dx%d (%.1f%% completed)',dimX,dimY,iter/nIter*100),
               asp=1,axes=FALSE,xlab=" ",ylab=" ") # 
          points(X[OutletPixel],Y[OutletPixel],pch=15,col=rnbw[catch[resort[OutletPixel]]])
          if (easyDraw == FALSE){
            for (i in AvailableNodesPlot){
              if (AA[i]<=(thrADraw)  & abs(X[i]-X[DownNode[i]])<=cellsize & abs(Y[i]-Y[DownNode[i]])<=cellsize ) {
                lines(c(X[i],X[DownNode[i]]),c(Y[i],Y[DownNode[i]]),lwd=0.5,col="#E0E0E0")}
            }
          }
          for (i in 1:Nnodes){
            if (!(i %in% OutletPixel)){
              if (AA[i]>(thrADraw)  & abs(X[i]-X[DownNode[i]])<=cellsize & abs(Y[i]-Y[DownNode[i]])<=cellsize ) {
                lines(c(X[i],X[DownNode[i]]),c(Y[i],Y[DownNode[i]]),lwd=0.5+4.5*(AA[i]/Nnodes/cellsize^2)^0.5,col=rnbw[catch[resort[i]]])}
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
  ######################
  ## EXPORT VARIABLES ##
  ######################
  
  #W <- as.dgCMatrix.spam(t(Wt)) # ensure compatibility with other functions
  W <- t(Wt)
  
  FD <- list(A=A*cellsize^2,W=W,downNode=DownNode,X=X,Y=Y,nNodes=Nnodes,outlet=OutletPixel,perm=pl)
  OCN <- list(FD=FD,dimX=dimX,dimY=dimY,cellsize=cellsize,nOutlet=nOutlet,periodicBoundaries=periodicBoundaries,
              expEnergy=expEnergy,coolingRate=coolingRate,typeInitialState=typeInitialState,nIter=nIter,initialNoCoolingPhase=initialNoCoolingPhase,
              energyInit=Energy_0)
  
  if (saveEnergy==TRUE) {OCN[["energy"]] <- Energy}
  if (saveExitFlag==TRUE) {OCN[["exitFlag"]] <- ExitFlag}
  if (saveN8==TRUE) {
    N8 <- list(W=W_N8)
    OCN[["N8"]] <- N8}
  if (saveN4==TRUE) {
    N4 <- list(W=W_N4)
    OCN[["N4"]] <- N4}
  
  invisible(OCN)
}


#########################################
## AUXILIARY FUNCTION initialstate_OCN ##
#########################################

initialstate_OCN <- function(dimX,dimY,nOutlet,outletSide,outletPos,typeInitialState){
  
  ## create initial state of the network
  flowDirStart <- matrix(data=0,nrow=dimY,ncol=dimX)
  
  if (typeInitialState=="H"){
    limit <- min(round(dimX/2),round(dimY/2))
    tmp <- matrix(c(1:limit,1:limit),limit,2)
    flowDirStart[tmp] <- 4
    tmp <- matrix(c(seq(limit,1,-1),(dimX-limit+1):dimX),limit,2)
    flowDirStart[tmp] <- 2
    tmp <- matrix(c((dimY-limit+1):dimY,seq(limit,1,-1)),limit,2)
    flowDirStart[tmp] <- 6
    tmp <- matrix(c((dimY-limit+1):dimY,(dimX-limit+1):dimX),limit,2)
    flowDirStart[tmp] <- 8
    for (i in (1:limit)){
      if ((dimX-i) >= (i+1)){
        flowDirStart[i,(i+1):(dimX-i)] <- 3
      }
    }
    for (i in (1:limit)){
      if ((dimX-i) >= (i+1)){
        flowDirStart[dimY+1-i,(i+1):(dimX-i)] <- 7
      }
    }
    for (i in (1:limit)){
      if ((dimY-i) >= (i+1)){
        flowDirStart[(i+1):(dimY-i),i] <- 5
      }
    }
    for (i in (1:limit)){
      if ((dimY-i) >= (i+1)){
        flowDirStart[(i+1):(dimY-i),dimX+1-i] <- 1
      }
    }
    #impose drains on the borders
    FirstSide <- outletSide[1]
    if (FirstSide=="N"){
      flowDirStart[1,(dimX-outletPos[1])] <- 5
      flowDirStart[1,(dimX-outletPos[1]):dimX] <- 1
      flowDirStart[,1] <- 7
      flowDirStart[,dimX] <- 7
      flowDirStart[dimY,1:outletPos[1]] <- 1
      flowDirStart[dimY,(outletPos[1]):dimX] <- 5
    } else if (FirstSide=="S"){
      flowDirStart[dimY,1:(dimX-outletPos[1])] <- 5
      flowDirStart[dimY,(dimX-outletPos[1]):dimX] <- 1
      flowDirStart[,1] <- 3
      flowDirStart[,dimX] <- 3
      flowDirStart[1,1:outletPos[1]] <- 1
      flowDirStart[1,(outletPos[1]):dimX] <- 5
    } else if (FirstSide=="W"){
      flowDirStart[1:(dimY-outletPos[1]),dimX] <- 3
      flowDirStart[(dimY-outletPos[1]):dimY,dimX] <- 7
      flowDirStart[1,] <- 5
      flowDirStart[dimY,] <- 5
      flowDirStart[1:outletPos[1],1] <- 7
      flowDirStart[(outletPos[1]):dimY,1] <- 3
    } else if (FirstSide=="E"){
      flowDirStart[1:(dimY-outletPos[1]),1] <- 3
      flowDirStart[(dimY-outletPos[1]):dimY,1] <- 7
      flowDirStart[1,] <- 1
      flowDirStart[dimY,] <- 1
      flowDirStart[1:outletPos[1],dimX] <- 7
      flowDirStart[(outletPos[1]):dimY,dimX] <- 3}
    
  }  else {
    
    if (typeInitialState=="I"){
      # flow direction matrix (defining a valley)
      FirstSide <-outletSide[1]
      if (FirstSide=="N"){
        if (outletPos[1]>1) {flowDirStart[,1:(outletPos[1]-1)] <- 1}
        if (dimX>outletPos[1]) {flowDirStart[,(outletPos[1]+1):dimX] <- 5}
        flowDirStart[,outletPos[1]] <- 7
      } else if (FirstSide=="E"){
        if (dimY>outletPos[1]) {flowDirStart[(outletPos[1]+1):dimY,] <- 3}
        if (outletPos[1]>1) {flowDirStart[1:(outletPos[1]-1),] <- 7}
        flowDirStart[outletPos[1],] <-1
      } else if (FirstSide=="S"){
        if (outletPos[1]>1) {flowDirStart[,1:(outletPos[1]-1)] <- 1}
        if (dimX>outletPos[1]) {flowDirStart[,(outletPos[1]+1):dimX] <- 5}
        flowDirStart[,outletPos[1]] <- 3
      } else if (FirstSide=="W"){
        if (dimY>outletPos[1]) {flowDirStart[(outletPos[1]+1):dimY,] <- 3}
        if (outletPos[1]>1) {flowDirStart[1:(outletPos[1]-1),] <- 7}
        flowDirStart[outletPos[1],] <- 5
      }
      # impose drains perpendicular to outlets
      # for (i in 1:nOutlet){
      #   if (outletSide[i]=="N"){
      #     flowDirStart[,outletPos[i]] <- 7 
      #   } else if (outletSide[i]=="E") {
      #     flowDirStart[outletPos[i],] <- 1 
      #   } else if (outletSide[i]=="S") {
      #     flowDirStart[,outletPos[i]] <- 3 
      #   } else if (outletSide[i]=="W") {
      #     flowDirStart[outletPos[i],] <- 5 
      #   }
      # }
      tmpX <- round(0.25*dimX); tmpY <- round(0.25*dimY)
      for (i in 1:nOutlet){
        if (outletSide[i]=="N"){
          lowX <- max(1,outletPos[i]-tmpX)
          uppX <- min(outletPos[i]+tmpX,dimX)
          flowDirStart[(dimY-tmpY+1):dimY,lowX:outletPos[i]] <- 1
          flowDirStart[(dimY-tmpY+1):dimY,outletPos[i]:uppX] <- 5
        } else if (outletSide[i]=="E") {
          lowY <- max(1,outletPos[i]-tmpY)
          uppY <- min(outletPos[i]+tmpY,dimY)
          flowDirStart[lowY:outletPos[i],(dimX-tmpX+1):dimX] <- 7
          flowDirStart[outletPos[i]:uppY,(dimX-tmpX+1):dimX] <- 3
        } else if (outletSide[i]=="S") {
          lowX <- max(1,outletPos[i]-tmpX)
          uppX <- min(outletPos[i]+tmpX,dimX)
          flowDirStart[1:tmpY,lowX:outletPos[i]] <- 1
          flowDirStart[1:tmpY,outletPos[i]:uppX] <- 5
        } else if (outletSide[i]=="W") {
          lowY <- max(1,outletPos[i]-tmpY)
          uppY <- min(outletPos[i]+tmpY,dimY)
          flowDirStart[lowY:outletPos[i],1:tmpX] <- 7
          flowDirStart[outletPos[i]:uppY,1:tmpX] <- 3
        }
      }
      for (i in 1:nOutlet){
        if (outletSide[i]=="N"){
          flowDirStart[(dimY-tmpY+1):dimY,outletPos[i]] <- 7 
        } else if (outletSide[i]=="E") {
          flowDirStart[outletPos[i],(dimX-tmpX+1):dimX] <- 1 
        } else if (outletSide[i]=="S") {
          flowDirStart[1:tmpY,outletPos[i]] <- 3 
        } else if (outletSide[i]=="W") {
          flowDirStart[outletPos[i],1:tmpX] <- 5 
        }
      } 
      
    } else if (typeInitialState=="T"){
      FirstSide <- outletSide[1]
      if (FirstSide=="N"){
        Transept <- round(0.5*dimY)
        if (outletPos[1]>1) {flowDirStart[(Transept+1):dimY,1:(outletPos[1]-1)] <- 1}
        if (dimX>outletPos[1]) {flowDirStart[(Transept+1):dimY,(outletPos[1]+1):dimX] <- 5}
        flowDirStart[(Transept+1):dimY,outletPos[1]] <- 7
        flowDirStart[1:Transept,] <- 7
      } else if (FirstSide=="E"){
        Transept <- round(0.5*dimX)
        if (dimY>outletPos[1]) {flowDirStart[(outletPos[1]+1):dimY,(Transept+1):dimX] <- 3}
        if (outletPos[1]>1) {flowDirStart[1:(outletPos[1]-1),(Transept+1):dimX] <- 7}
        flowDirStart[outletPos[1],(Transept+1):dimX] <- 1
        flowDirStart[,1:Transept] <- 1
      } else if (FirstSide=="S"){
        Transept <- round(0.5*dimY)
        if (outletPos[1]>1) {flowDirStart[1:Transept,1:(outletPos[1]-1)] <- 1}
        if (dimX>outletPos[1]) {flowDirStart[1:Transept,(outletPos[1]+1):dimX] <- 5}
        flowDirStart[1:Transept,outletPos[1]] <- 3
        flowDirStart[(Transept+1):dimY,] <- 3
      } else if (FirstSide=="W"){
        Transept <- round(0.5*dimX)
        if (dimY>outletPos[1]) {flowDirStart[(outletPos[1]+1):dimY,1:Transept] <- 3}
        if (outletPos[1]>1) {flowDirStart[1:(outletPos[1]-1),1:Transept] <- 7}
        flowDirStart[outletPos[1],1:Transept] <- 5
        flowDirStart[,(Transept+1):dimX] <- 5
      }
      # impose drains perpendicular to outlets
      tmpX <- round(0.25*dimX); tmpY <- round(0.25*dimY)
      for (i in 1:nOutlet){
        if (outletSide[i]=="N"){
          lowX <- max(1,outletPos[i]-tmpX)
          uppX <- min(outletPos[i]+tmpX,dimX)
          flowDirStart[(dimY-tmpY+1):dimY,lowX:outletPos[i]] <- 1
          flowDirStart[(dimY-tmpY+1):dimY,outletPos[i]:uppX] <- 5
        } else if (outletSide[i]=="E") {
          lowY <- max(1,outletPos[i]-tmpY)
          uppY <- min(outletPos[i]+tmpY,dimY)
          flowDirStart[lowY:outletPos[i],(dimX-tmpX+1):dimX] <- 7
          flowDirStart[outletPos[i]:uppY,(dimX-tmpX+1):dimX] <- 3
        } else if (outletSide[i]=="S") {
          lowX <- max(1,outletPos[i]-tmpX)
          uppX <- min(outletPos[i]+tmpX,dimX)
          flowDirStart[1:tmpY,lowX:outletPos[i]] <- 1
          flowDirStart[1:tmpY,outletPos[i]:uppX] <- 5
        } else if (outletSide[i]=="W") {
          lowY <- max(1,outletPos[i]-tmpY)
          uppY <- min(outletPos[i]+tmpY,dimY)
          flowDirStart[lowY:outletPos[i],1:tmpX] <- 7
          flowDirStart[outletPos[i]:uppY,1:tmpX] <- 3
        }
      }
      for (i in 1:nOutlet){
        if (outletSide[i]=="N"){
          flowDirStart[(dimY-tmpY+1):dimY,outletPos[i]] <- 7 
        } else if (outletSide[i]=="E") {
          flowDirStart[outletPos[i],(dimX-tmpX+1):dimX] <- 1 
        } else if (outletSide[i]=="S") {
          flowDirStart[1:tmpY,outletPos[i]] <- 3 
        } else if (outletSide[i]=="W") {
          flowDirStart[outletPos[i],1:tmpX] <- 5 
        }
      }      
    } else if (typeInitialState=="V") {
      FirstSide <- outletSide[1]
      if (FirstSide=="N"){
        for (i in 1:min(outletPos[1],dimY)) {flowDirStart[dimY+1-i,outletPos[1]+1-i]=8}
        if (dimX>outletPos[1]) {for (i in 1:min(dimX-outletPos[1],dimY-1)) {flowDirStart[dimY-i,outletPos[1]+i]=6}}
        if (outletPos[1]>1) {for (i in 1:min(outletPos[1]-1,dimY)) {flowDirStart[dimY+1-i,1:(outletPos[1]-i)]=1}}
        if (dimX>outletPos[1]) {for (i in 1:min(dimY,dimX-outletPos[1])) {flowDirStart[dimY+1-i,(outletPos[1]+i):dimX]=5}}
        flowDirStart[flowDirStart==0]=7
      }
      if (FirstSide=="E"){
        for (i in 1:min(outletPos[1],dimY)) {flowDirStart[i,dimX-outletPos[1]+i]=8}
        if (dimY>outletPos[1]) {for (i in 1:min(dimY-outletPos[1],dimX-1)) {flowDirStart[outletPos[1]+i,dimX-i]=2}}
        if (outletPos[1]>1) {for (i in 1:min((outletPos[1]-1),dimY)) {flowDirStart[i,(dimX-outletPos[1]+i+1):dimX]=7}}
        if (dimY>outletPos[1]) {for (i in 1:min(dimY-outletPos[1],dimX)) {flowDirStart[(outletPos[1]+i):dimY,dimX+1-i]=3}}
        flowDirStart[flowDirStart==0]=1
      }
      if (FirstSide=="S"){
        for (i in 1:min(outletPos[1],dimY)) {flowDirStart[i,outletPos[1]+1-i]=2}
        if (dimX>outletPos[1]) {for (i in 1:min(dimX-outletPos[1],dimY-1)) {flowDirStart[i+1,outletPos[1]+i]=4}}
        if (outletPos[1]>1) {for (i in 1:min(outletPos[1]-1,dimY)) {flowDirStart[i,1:(outletPos[1]-i)]=1}}
        if (dimX>outletPos[1]) {for (i in 1:min(dimX-outletPos[1],dimY)) {flowDirStart[i,(outletPos[1]+i):dimX]=5}}
        flowDirStart[flowDirStart==0]=3
      }
      if (FirstSide=="W"){
        for (i in 1:min(outletPos[1],dimY)) {flowDirStart[i,outletPos[1]-i+1]=6}
        if (dimY>outletPos[1]) {for (i in 1:min(dimY-outletPos[1],dimX-1)) {flowDirStart[outletPos[1]+i,i+1]=4}}
        if (outletPos[1]>1) {for (i in 1:min(outletPos[1]-1,dimY)) {flowDirStart[i,1:(outletPos[1]-i)]=7}}
        if (dimY>outletPos[1]) {for (i in 1:min(dimY-outletPos[1],dimX)) {flowDirStart[(outletPos[1]+i):dimY,i]=3}}
        flowDirStart[flowDirStart==0]=5
      }
      
      # impose drains with V shape
      # for (i in 1:nOutlet){
      #   if (outletSide[i]=="N"){
      #     for (j in 1:min(outletPos[i],dimY)) {flowDirStart[dimY+1-j,outletPos[i]+1-j]=8}
      #     if (dimX>outletPos[i]) {for (j in 1:min(dimX-outletPos[i],dimY-1)) {flowDirStart[dimY-j,outletPos[i]+j]=6}}
      #   } else if (outletSide[i]=="E") {
      #     for (j in 1:outletPos[i]) {flowDirStart[j,dimX-outletPos[i]+j]=8}
      #     if (dimY>outletPos[i]) {for (j in 1:min(dimY-outletPos[i],dimX-1)) {flowDirStart[outletPos[i]+j,dimX-j]=2}}
      #   } else if (outletSide[i]=="S") {
      #     for (j in 1:min(outletPos[i],dimY)) {flowDirStart[j,outletPos[i]+1-j]=2}
      #     if (dimX>outletPos[i]) {for (j in 1:min(dimX-outletPos[i],dimY-1)) {flowDirStart[j+1,outletPos[i]+j]=4}}
      #   } else if (outletSide[i]=="W") {
      #     for (j in 1:outletPos[i]) {flowDirStart[j,outletPos[i]-j+1]=6}
      #     if (dimY>outletPos[i]) {for (j in 1:min(dimY-outletPos[i],dimX-1)) {flowDirStart[outletPos[i]+j,j+1]=4}}
      #   }
      # }
      if (nOutlet > 1){ ## BUG WHEN OUTLET IS ON LONG SIDE OF A RECTANGULAR LATTICE
        tmpX <- round(0.25*dimX); tmpY <- round(0.25*dimY)
        for (i in 1:nOutlet){
          if (outletSide[i]=="N"){
            lowX <- max(1,outletPos[i]-tmpX)
            uppX <- min(outletPos[i]+tmpX,dimX)
            flowDirStart[(dimY-tmpY):dimY,lowX:uppX] <- 0
            if (outletPos[i]>1) {for (j in 1:(outletPos[i]-lowX)) {flowDirStart[dimY+1-j,lowX:(outletPos[i]-j)]=1}}
            if (dimX>outletPos[i]) {for (j in 1:(uppX-outletPos[i])) {flowDirStart[dimY+1-j,(outletPos[i]+j):uppX]=5}}
            flowDirStart[flowDirStart==0]=7
            
          } else if (outletSide[i]=="S"){
            lowX <- max(1,outletPos[i]-tmpX)
            uppX <- min(outletPos[i]+tmpX,dimX)
            flowDirStart[1:(tmpY+1),lowX:uppX] <- 0
            if (outletPos[i]>1) {for (j in 1:(outletPos[i]-lowX)) {flowDirStart[j,lowX:(outletPos[i]-j)]=1}}
            if (dimX>outletPos[i]) {for (j in 1:(uppX-outletPos[i])) {flowDirStart[j,(outletPos[i]+j):uppX]=5}}
            flowDirStart[flowDirStart==0]=3
            
          } else if (outletSide[i]=="E"){
            lowY <- max(1,outletPos[i]-tmpY)
            uppY <- min(outletPos[i]+tmpY,dimY)
            flowDirStart[lowY:uppY,(dimX-tmpX):dimX] <- 0
            if (outletPos[i]>1) {for (j in 1:(outletPos[i]-lowY)) {flowDirStart[lowY-1+j,(dimX-outletPos[i]+lowY+j-1):dimX]=7}}
            if (dimY>outletPos[i]) {for (j in 1:(uppY-outletPos[i])) {flowDirStart[outletPos[i]+j,(dimX+1-j):dimX]=3}}
            flowDirStart[flowDirStart==0]=1
          }  else if (outletSide[i]=="W"){
            lowY <- max(1,outletPos[i]-tmpY)
            uppY <- min(outletPos[i]+tmpY,dimY)
            flowDirStart[lowY:uppY,1:(tmpX+1)] <- 0
            if (outletPos[i]>1) {for (j in 1:(outletPos[i]-lowY)) {flowDirStart[lowY-1+j,1:(outletPos[i]-lowY-j+1)]=7}}
            if (dimY>outletPos[i]) {for (j in 1:(uppY-outletPos[i])) {flowDirStart[outletPos[i]+j,1:j]=3}}
            flowDirStart[flowDirStart==0]=5
          }
        }
        for (i in 1:nOutlet){
          if (outletSide[i]=="N"){
            lowX <- max(1,outletPos[i]-tmpX)
            uppX <- min(outletPos[i]+tmpX,dimX)
            if (outletPos[i]>1) {for (j in 1:(outletPos[i]-lowX+1)) {flowDirStart[dimY+1-j,outletPos[i]+1-j]=8}}
            if (dimX>outletPos[i]) {for (j in 1:(uppX-outletPos[i])) {flowDirStart[dimY-j,outletPos[i]+j]=6}}
            
          } else if (outletSide[i]=="S"){
            lowX <- max(1,outletPos[i]-tmpX)
            uppX <- min(outletPos[i]+tmpX,dimX)
            if (outletPos[i]>1) {for (j in 1:(outletPos[i]-lowX+1)) {flowDirStart[j,outletPos[i]+1-j]=2}}
            if (dimX>outletPos[i]) {for (j in 1:(uppX-outletPos[i])) {flowDirStart[j+1,outletPos[i]+j]=4}}
            
          } else if (outletSide[i]=="E"){
            lowY <- max(1,outletPos[i]-tmpY)
            uppY <- min(outletPos[i]+tmpY,dimY)
            if (outletPos[i]>1) {for (j in 1:(outletPos[i]-lowY)) {flowDirStart[lowY-1+j,dimX-outletPos[i]+lowY+j-1]=8}}
            if (dimY>outletPos[i]) {for (j in 1:(uppY-outletPos[i])) {flowDirStart[outletPos[i]+j,dimX-j]=2}}
          }  else if (outletSide[i]=="W"){
            lowY <- max(1,outletPos[i]-tmpY)
            uppY <- min(outletPos[i]+tmpY,dimY)
            if (outletPos[i]>1) {for (j in 1:(outletPos[i]-lowY)) {flowDirStart[lowY-1+j,outletPos[i]-lowY-j+2]=6}}
            if (dimY>outletPos[i]) {for (j in 1:(uppY-outletPos[i])) {flowDirStart[outletPos[i]+j,j+1]=4}}
          }
        }
      }
    } else {
      stop("Invalid initial state")
    }
  }
  invisible(flowDirStart)
}


initial_permutation <- function(DownNode){
  
  Outlet <- which(DownNode==0)
  NodesToExplore <- Outlet # start from outlets
  reverse_perm <- numeric(length(DownNode)) # build permutation vector from outlets to headwaters, then flip it
  
  k <- 0
  while (length(NodesToExplore)>0){ # continue until all the network has been explored
    k <- k + 1
    node <- NodesToExplore[1] # explore a node
    reverse_perm[k] <- node # assign position in the permutation vector
    NodesToExplore <- NodesToExplore[-1] # remove explored node
    UpNodes <- which(DownNode==node) # find nodes upstream of node
    while (length(UpNodes)>0){ # continue upstream until a headwater is found
      k <- k + 1
      node <- UpNodes[1] # explore first upstream node
      reverse_perm[k] <- node
      if (length(UpNodes)>1){ # if there is a bifurcation upstream, add the other upstream connections at the top of NodesToExplore
        NodesToExplore <- c(UpNodes[2:length(UpNodes)],NodesToExplore)
      }
      UpNodes <- which(DownNode==node)
    }
  }
  
  perm <- reverse_perm[length(DownNode):1] # flip permutation
  
  OutList = list(perm=perm,noDAG=0)
  
  invisible(OutList)
}


