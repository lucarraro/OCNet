
create_OCN <- function(dimX,dimY,
                       N_outlet=1,
                       OutletSide="S", # vector of sides where outlets are located
                       OutletPos=round(dimX/3),
                       PeriodicBoundaries=FALSE, # if TRUE, reflecting boundaries are applied. If FALSE, boundaries are non-reflecting.
                       type_initialstate=NULL,
                       FlowDirStart=NULL,
                       ExpEnergy=0.5, # energy \propto Q*deltaH, Q propto A, deltaH propto A^-0.5
                       cellsize=1,
                       xllcorner=0.5*cellsize,
                       yllcorner=0.5*cellsize,
                       N_iter=30*dimX*dimY,
                       N_updates=50, # number of times an update is shown
                       InitialNoCoolingPhase=0,
                       CoolingRate=50,
                       ShowIntermediatePlots=FALSE,
                       A_thr_draw=0.002*dimX*dimY*cellsize^2,
                       SaveEnergy=FALSE,
                       SaveExitFlag=FALSE,
                       SaveN8=FALSE,
                       SaveN4=FALSE,
                       DisplayUpdates=1){
 
  
  if (dimX<2 | dimY<2) stop("Dimensions too small.")
  if (dimX*dimY > 1000) {
    dim <- sqrt(dimX*dimY)
    estTime <- (7.988e-1 - 1.266e-1*dim + 7.198e-3*dim^2 - 6.807e-5*dim^3 +  1.372e-6*dim^4) / (30*dimX*dimY) * N_iter
    unitTime <- "seconds"
    if (estTime >= 60 && estTime < 3600 ){
      estTime <- estTime/60
      unitTime <- "minutes"
    } else if (estTime >= 3600) {
      estTime <- estTime/3600
    unitTime <- "hours"}
    
    if (DisplayUpdates != 0){
      cat("create_OCN is running...\n")
      cat(sprintf("Estimated duration: %.2f %s \n",estTime,unitTime))
      cat("Note that the above estimate is only based on the choice of parameters dimX, dimY and N_iter, and not on processor performance.\n")
      cat("\n")
    }
  }
  
  t0 <- Sys.time()
  if (DisplayUpdates==2){cat('Initializing...\n')}
  
  ################################## 
  ## DEFINE INITIAL NETWORK STATE ##  
  ##################################  
  
  if ((N_outlet=="All")==TRUE && is.null(type_initialstate)==TRUE){
    type_initialstate <- "H"
  } else {
    if (is.null(type_initialstate)==TRUE){
      type_initialstate <- "I"
    }
  }
  
  # set OutletSide, OutletPos in accordance with the input given 
  if ((N_outlet=="All")==TRUE){
    OutletSide <- c(rep("S",dimX),rep("E",dimY-2),rep("N",dimX),rep("W",dimY-2))
    OutletPos <- c(1:dimX,2:(dimY-1),seq(dimX,1,-1),seq((dimY-1),2,-1))
    N_outlet <- 2*(dimX+dimY-2)
  } else if (N_outlet>1 && length(OutletSide)==1 && length(OutletPos)==1 && all(OutletSide=="S") && all(OutletPos==round(dimX/3))) {
    BorderPixels <- 1:(2*(dimX+dimY-2))
    AllSides <- c(rep("S",dimX),rep("E",dimY-2),rep("N",dimX),rep("W",dimY-2))
    AllPos <- c(1:dimX,2:(dimY-1),1:dimX,2:(dimY-1))
    outlet_ID <- sample(BorderPixels,N_outlet)
    OutletSide <- AllSides[outlet_ID]
    OutletPos <- AllPos[outlet_ID]
    # for (i in 1:N_outlet){
    #   OutletSide[i] <- sample(AllSides,1)
    #   if (OutletSide[i]=="N" || OutletSide[i]=="S"){
    #     OutletPos[i] <- sample(1:dimX,1) 
    #   } else {OutletPos[i] <- sample(1:dimY,1)}
    # }
  }
  
  # check whether input is correct
  if(N_outlet != length(OutletSide) || N_outlet != length(OutletPos) ){
    stop('Length of OutletSide and/or OutletPos is inconsistent with N_outlet')}
  
  if (!all(OutletSide %in% c("N","W","S","E"))){
    stop('Invalid OutletSide')}
  
  for (i in 1:N_outlet){
    if (OutletSide[i] %in% c("N","S")){
      if (!(OutletPos[i] %in% 1:dimX)){
        stop('Invalid OutletPos')} 
    } else if (!(OutletPos[i] %in% 1:dimY)){
      stop('Invalid OutletPos')}
  }
  
  if (!(type_initialstate %in% c("I","T","V","H","R"))){
    stop('Invalid type_initialstate')}
  
  
  # define domain coordinates
  X <- rep(seq(xllcorner,xllcorner+(dimX-1)*cellsize,cellsize), each = dimY)  
  Y <- rep(seq(yllcorner,yllcorner+(dimY-1)*cellsize,cellsize),dimX)
  
  # find list of possible neighbouring pixels
  Nnodes <- length(X)
  movement <- matrix(c(0,-1,-1,-1,0,1,1,1,1,1,0,-1,-1,-1,0,1),nrow=2,byrow=TRUE)
  NeighbouringNodes <- vector("list", Nnodes)
  if (SaveN8) {W_N8 <- spam(0, Nnodes,Nnodes)}  # null matrices
  if (SaveN4) {W_N4 <- spam(0, Nnodes,Nnodes)}  
  cont_node <- 0
  for (cc in 1:dimX) {
    for (rr in 1:dimY) {
      cont_node <- cont_node + 1
      neigh_r <- rr+movement[1,]
      neigh_c <- cc+movement[2,]
      if (PeriodicBoundaries==TRUE){
        neigh_r[neigh_r==0] <- dimY
        neigh_c[neigh_c==0] <- dimX
        neigh_r[neigh_r>dimY] <- 1
        neigh_c[neigh_c>dimX] <- 1
      }
      NotAboundary <- neigh_r>0 & neigh_r<=dimY & neigh_c>0 & neigh_c<=dimX # only effective when PeriodicBoundaries=FALSE
      NotAboundary_N4 <- neigh_r>0 & neigh_r<=dimY & neigh_c>0 & neigh_c<=dimX & (((1:8) %% 2)==1)
      NeighbouringNodes[[cont_node]] <- neigh_r[NotAboundary] + (neigh_c[NotAboundary]-1)*dimY
      if (SaveN8){W_N8[cont_node,neigh_r[NotAboundary] + (neigh_c[NotAboundary]-1)*dimY]=1}
      if (SaveN4){W_N4[cont_node,neigh_r[NotAboundary_N4] + (neigh_c[NotAboundary_N4]-1)*dimY]=1}
    }
  } 
  
  # find outlet position
  Outlet_row <- numeric(N_outlet)
  Outlet_col <- numeric(N_outlet)
  for (i in 1:N_outlet){
    if (OutletSide[i]=="N"){
      Outlet_row[i] <- dimY
      Outlet_col[i] <- OutletPos[i]
    } else if (OutletSide[i]=="E") {
      Outlet_row[i] <- OutletPos[i]
      Outlet_col[i] <- dimX
    } else if (OutletSide[i]=="S") {
      Outlet_row[i] <- 1
      Outlet_col[i] <- OutletPos[i]
    } else if (OutletSide[i]=="W") {
      Outlet_row[i] <- OutletPos[i]
      Outlet_col[i] <- 1
    }
  }
  
  # find flow direction matrix given initial network state
  if (is.null(FlowDirStart)==TRUE){
    FlowDirStart <- initialstate_OCN(dimX,dimY,N_outlet,OutletSide,OutletPos,type_initialstate)}
  
  # attribute flow direction =0 to outlet pixels
  for (i in 1:N_outlet){
    if (OutletSide[i]=="N"){
      FlowDirStart[dimY,OutletPos[i]] <- 0 
      Outlet_row[i] <- dimY
      Outlet_col[i] <- OutletPos[i]
    } else if (OutletSide[i]=="E") {
      FlowDirStart[OutletPos[i],dimX] <- 0 
      Outlet_row[i] <- OutletPos[i]
      Outlet_col[i] <- dimX
    } else if (OutletSide[i]=="S") {
      FlowDirStart[1,OutletPos[i]] <- 0
      Outlet_row[i] <- 1
      Outlet_col[i] <- OutletPos[i]
    } else if (OutletSide[i]=="W") {
      FlowDirStart[OutletPos[i],1] <- 0
      Outlet_row[i] <- OutletPos[i]
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
      dir <- FlowDirStart[rr,cc]
      if (dir>0) {
        DownPixel <- c(rr,cc)+movement[,dir]
        # if a custom FlowDirStart is provided, allow flow to cross boundaries
        if (PeriodicBoundaries==TRUE){
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
  AvailableNodes <- setdiff(1:Nnodes,OutletPixel)
  
  pl <- initial_permutation(DownNode)  # calculate permutation vector
  
  pas <- permuteAddSolve(Wt, pl$perm, 1L, numeric(Nnodes), ExpEnergy)
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
  Energy <- numeric(N_iter)
  Energy_0 <- pas[[2]]
  Energy[1] <-Energy_0
  Temperature <- c(Energy[1]+numeric(InitialNoCoolingPhase*N_iter),Energy[1]*exp(-CoolingRate*(1:(N_iter-InitialNoCoolingPhase*N_iter))/N_iter))
  savetime <- numeric(N_iter)
  ExitFlag <- numeric(N_iter)
  
  # plot initial state
  if (ShowIntermediatePlots==TRUE){
    AA <- A*cellsize^2
    if (N_outlet > 1){
      rnbw <- hcl.colors(N_outlet,palette="Dark 3")
      rnbw <- c(rnbw[(round(N_outlet/2)+1):N_outlet],rnbw[1:(round(N_outlet/2)+1)])
    } else {rnbw <- hcl.colors(3,palette="Dark 3")
    rnbw <- rnbw[3]}
    resort <- invPerm(pl$perm)
    catch <- numeric(Nnodes)
    for (o in 1:N_outlet){
      catch[(resort[OutletPixel[o]]-AA[OutletPixel[o]]/cellsize^2+1):resort[OutletPixel[o]]] <- o
    }
    par(bty="n")
    plot(c(min(X),max(X)),c(min(Y),max(Y)),main=sprintf('OCN %dx%d (initial state)',dimX,dimY),
         type="n",asp=1,axes=FALSE,xlab="",ylab="") # 
    points(X[OutletPixel],Y[OutletPixel],pch=15,col=rnbw[catch[resort[OutletPixel]]])
    for (i in AvailableNodes){
      if (AA[i]<=A_thr_draw & abs(X[i]-X[DownNode[i]])<=cellsize & abs(Y[i]-Y[DownNode[i]])<=cellsize ) {
        lines(c(X[i],X[DownNode[i]]),c(Y[i],Y[DownNode[i]]),lwd=0.5,col="#E0E0E0")}
    }
    for (i in 1:Nnodes){
      if (!(i %in% OutletPixel)){
      if (AA[i]>A_thr_draw & abs(X[i]-X[DownNode[i]])<=cellsize & abs(Y[i]-Y[DownNode[i]])<=cellsize ) {
        lines(c(X[i],X[DownNode[i]]),c(Y[i],Y[DownNode[i]]),lwd=0.5+4.5*(AA[i]/Nnodes/cellsize^2)^0.5,col=rnbw[catch[resort[i]]])}
    }}
  }
  
  #Rprof("\\\\eawag/userdata/carrarlu/Desktop/Luca/OCN/R/spam_stuff/Rprof100.out")
  if (DisplayUpdates == 2){
    cat(sprintf('Initialization completed. Elapsed time is %.2f s \n',difftime(Sys.time(),t0,units='secs'))) 
    cat('Search algorithm has started.\n')
    cat('\n')}
  
  pl <- as.integer(pl$perm)
  flag <- 0
  
  # simulated annealing algorithm
  for (iter in 2:N_iter) {
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
      
        pas <- allinoneF( as.integer(N_outlet), pl, Wt, DownNode, node,
                        down_new, A[node], ExpEnergy)
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
    if (DisplayUpdates==2){
      if (iter %% round(N_iter/N_updates)==0){
        cat(sprintf('%.1f%% completed - Elapsed time: %.2f s - %11s - Energy; %0.f \n',
                    iter/N_iter*100,difftime(Sys.time(),t0,units='secs'),format(Sys.time(),"%b%d %H:%M"),Energy[iter]))
        }}
    # plot update
    #if (ShowIntermediatePlots==TRUE){
    #par(bty="n")
    
    if (iter %% round(N_iter/N_updates)==0){
      if (ShowIntermediatePlots==TRUE){ 
        AA <- A*cellsize^2 
        resort <- invPerm(pl)
        catch <- numeric(Nnodes)
        for (o in 1:N_outlet){
          catch[(resort[OutletPixel[o]]-AA[OutletPixel[o]]/cellsize^2+1):resort[OutletPixel[o]]] <- o
        }
        plot(c(min(X),max(X)),c(min(Y),max(Y)),type="n",main=sprintf('OCN %dx%d (%.1f%% completed)',dimX,dimY,iter/N_iter*100),
             asp=1,axes=FALSE,xlab=" ",ylab=" ") # 
        points(X[OutletPixel],Y[OutletPixel],pch=15,col=rnbw[catch[resort[OutletPixel]]])
        for (i in AvailableNodes){
          if (AA[i]<=(A_thr_draw)  & abs(X[i]-X[DownNode[i]])<=cellsize & abs(Y[i]-Y[DownNode[i]])<=cellsize ) {
            lines(c(X[i],X[DownNode[i]]),c(Y[i],Y[DownNode[i]]),lwd=0.5,col="#E0E0E0")}
        }
        for (i in 1:Nnodes){
          if (!(i %in% OutletPixel)){
          if (AA[i]>(A_thr_draw)  & abs(X[i]-X[DownNode[i]])<=cellsize & abs(Y[i]-Y[DownNode[i]])<=cellsize ) {
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
  
  ######################
  ## EXPORT VARIABLES ##
  ######################
  
  #W <- as.dgCMatrix.spam(t(Wt)) # ensure compatibility with other functions
  W <- t(Wt)
  
  FD <- list(A=A*cellsize^2,W=W,DownNode=DownNode,X=X,Y=Y,Nnodes=Nnodes,Outlet=OutletPixel)
  OCN <- list(FD=FD,dimX=dimX,dimY=dimY,cellsize=cellsize,N_outlet=N_outlet,PeriodicBoundaries=PeriodicBoundaries,
              ExpEnergy=ExpEnergy,CoolingRate=CoolingRate,type_initialstate=type_initialstate,N_iter=N_iter,InitialNoCoolingPhase=InitialNoCoolingPhase)
  
  if (SaveEnergy==TRUE) {OCN[["Energy"]] <- Energy}
  if (SaveExitFlag==TRUE) {OCN[["ExitFlag"]] <- ExitFlag}
  if (SaveN8==TRUE) {
    N8 <- list(W=W_N8)
    OCN[["N8"]] <- N8}
  if (SaveN4==TRUE) {
    N4 <- list(W=W_N4)
    OCN[["N4"]] <- N4}
  
  return(OCN)
}


#########################################
## AUXILIARY FUNCTION initialstate_OCN ##
#########################################

initialstate_OCN <- function(dimX,dimY,N_outlet,OutletSide,OutletPos,type_initialstate){

  ## create initial state of the network
  FlowDirStart <- matrix(data=0,nrow=dimY,ncol=dimX)
  
  if (type_initialstate=="H"){
    limit <- min(round(dimX/2),round(dimY/2))
    tmp <- matrix(c(1:limit,1:limit),limit,2)
    FlowDirStart[tmp] <- 4
    tmp <- matrix(c(seq(limit,1,-1),(dimX-limit+1):dimX),limit,2)
    FlowDirStart[tmp] <- 2
    tmp <- matrix(c((dimY-limit+1):dimY,seq(limit,1,-1)),limit,2)
    FlowDirStart[tmp] <- 6
    tmp <- matrix(c((dimY-limit+1):dimY,(dimX-limit+1):dimX),limit,2)
    FlowDirStart[tmp] <- 8
    for (i in (1:limit)){
      if ((dimX-i) >= (i+1)){
        FlowDirStart[i,(i+1):(dimX-i)] <- 3
      }
    }
    for (i in (1:limit)){
      if ((dimX-i) >= (i+1)){
        FlowDirStart[dimY+1-i,(i+1):(dimX-i)] <- 7
      }
    }
    for (i in (1:limit)){
      if ((dimY-i) >= (i+1)){
        FlowDirStart[(i+1):(dimY-i),i] <- 5
      }
    }
    for (i in (1:limit)){
      if ((dimY-i) >= (i+1)){
        FlowDirStart[(i+1):(dimY-i),dimX+1-i] <- 1
      }
    }
    #impose drains on the borders
    FirstSide <- OutletSide[1]
    if (FirstSide=="N"){
      FlowDirStart[1,(dimX-OutletPos[1])] <- 5
      FlowDirStart[1,(dimX-OutletPos[1]):dimX] <- 1
      FlowDirStart[,1] <- 7
      FlowDirStart[,dimX] <- 7
      FlowDirStart[dimY,1:OutletPos[1]] <- 1
      FlowDirStart[dimY,(OutletPos[1]):dimX] <- 5
    } else if (FirstSide=="S"){
      FlowDirStart[dimY,1:(dimX-OutletPos[1])] <- 5
      FlowDirStart[dimY,(dimX-OutletPos[1]):dimX] <- 1
      FlowDirStart[,1] <- 3
      FlowDirStart[,dimX] <- 3
      FlowDirStart[1,1:OutletPos[1]] <- 1
      FlowDirStart[1,(OutletPos[1]):dimX] <- 5
    } else if (FirstSide=="W"){
      FlowDirStart[1:(dimY-OutletPos[1]),dimX] <- 3
      FlowDirStart[(dimY-OutletPos[1]):dimY,dimX] <- 7
      FlowDirStart[1,] <- 5
      FlowDirStart[dimY,] <- 5
      FlowDirStart[1:OutletPos[1],1] <- 7
      FlowDirStart[(OutletPos[1]):dimY,1] <- 3
    } else if (FirstSide=="E"){
      FlowDirStart[1:(dimY-OutletPos[1]),1] <- 3
      FlowDirStart[(dimY-OutletPos[1]):dimY,1] <- 7
      FlowDirStart[1,] <- 1
      FlowDirStart[dimY,] <- 1
      FlowDirStart[1:OutletPos[1],dimX] <- 7
      FlowDirStart[(OutletPos[1]):dimY,dimX] <- 3}
    
  }  else {
    
    if (type_initialstate=="I"){
      # flow direction matrix (defining a valley)
      FirstSide <-OutletSide[1]
      if (FirstSide=="N"){
        if (OutletPos[1]>1) {FlowDirStart[,1:(OutletPos[1]-1)] <- 1}
        if (dimX>OutletPos[1]) {FlowDirStart[,(OutletPos[1]+1):dimX] <- 5}
        FlowDirStart[,OutletPos[1]] <- 7
      } else if (FirstSide=="E"){
        if (dimY>OutletPos[1]) {FlowDirStart[(OutletPos[1]+1):dimY,] <- 3}
        if (OutletPos[1]>1) {FlowDirStart[1:(OutletPos[1]-1),] <- 7}
        FlowDirStart[OutletPos[1],] <-1
      } else if (FirstSide=="S"){
        if (OutletPos[1]>1) {FlowDirStart[,1:(OutletPos[1]-1)] <- 1}
        if (dimX>OutletPos[1]) {FlowDirStart[,(OutletPos[1]+1):dimX] <- 5}
        FlowDirStart[,OutletPos[1]] <- 3
      } else if (FirstSide=="W"){
        if (dimY>OutletPos[1]) {FlowDirStart[(OutletPos[1]+1):dimY,] <- 3}
        if (OutletPos[1]>1) {FlowDirStart[1:(OutletPos[1]-1),] <- 7}
        FlowDirStart[OutletPos[1],] <- 5
      }
      # impose drains perpendicular to outlets
      # for (i in 1:N_outlet){
      #   if (OutletSide[i]=="N"){
      #     FlowDirStart[,OutletPos[i]] <- 7 
      #   } else if (OutletSide[i]=="E") {
      #     FlowDirStart[OutletPos[i],] <- 1 
      #   } else if (OutletSide[i]=="S") {
      #     FlowDirStart[,OutletPos[i]] <- 3 
      #   } else if (OutletSide[i]=="W") {
      #     FlowDirStart[OutletPos[i],] <- 5 
      #   }
      # }
      tmpX <- round(0.15*dimX); tmpY <- round(0.15*dimY)
      for (i in 1:N_outlet){
         if (OutletSide[i]=="N"){
           lowX <- max(1,OutletPos[i]-tmpX)
           uppX <- min(OutletPos[i]+tmpX,dimX)
           FlowDirStart[(dimY-tmpY+1):dimY,lowX:OutletPos[i]] <- 1
           FlowDirStart[(dimY-tmpY+1):dimY,OutletPos[i]:uppX] <- 5
         } else if (OutletSide[i]=="E") {
           lowY <- max(1,OutletPos[i]-tmpY)
           uppY <- min(OutletPos[i]+tmpY,dimY)
           FlowDirStart[lowY:OutletPos[i],(dimX-tmpX+1):dimX] <- 7
           FlowDirStart[OutletPos[i]:uppY,(dimX-tmpX+1):dimX] <- 3
         } else if (OutletSide[i]=="S") {
           lowX <- max(1,OutletPos[i]-tmpX)
           uppX <- min(OutletPos[i]+tmpX,dimX)
           FlowDirStart[1:tmpY,lowX:OutletPos[i]] <- 1
           FlowDirStart[1:tmpY,OutletPos[i]:uppX] <- 5
         } else if (OutletSide[i]=="W") {
           lowY <- max(1,OutletPos[i]-tmpY)
           uppY <- min(OutletPos[i]+tmpY,dimY)
           FlowDirStart[lowY:OutletPos[i],1:tmpX] <- 7
           FlowDirStart[OutletPos[i]:uppY,1:tmpX] <- 3
         }
       }
      for (i in 1:N_outlet){
        if (OutletSide[i]=="N"){
          FlowDirStart[(dimY-tmpY+1):dimY,OutletPos[i]] <- 7 
        } else if (OutletSide[i]=="E") {
          FlowDirStart[OutletPos[i],(dimX-tmpX+1):dimX] <- 1 
        } else if (OutletSide[i]=="S") {
          FlowDirStart[1:tmpY,OutletPos[i]] <- 3 
        } else if (OutletSide[i]=="W") {
          FlowDirStart[OutletPos[i],1:tmpX] <- 5 
        }
      } 
      
    } else if (type_initialstate=="T"){
      FirstSide <- OutletSide[1]
      if (FirstSide=="N"){
        Transept <- round(0.5*dimY)
        if (OutletPos[1]>1) {FlowDirStart[(Transept+1):dimY,1:(OutletPos[1]-1)] <- 1}
        if (dimX>OutletPos[1]) {FlowDirStart[(Transept+1):dimY,(OutletPos[1]+1):dimX] <- 5}
        FlowDirStart[(Transept+1):dimY,OutletPos[1]] <- 7
        FlowDirStart[1:Transept,] <- 7
      } else if (FirstSide=="E"){
        Transept <- round(0.5*dimX)
        if (dimY>OutletPos[1]) {FlowDirStart[(OutletPos[1]+1):dimY,(Transept+1):dimX] <- 3}
        if (OutletPos[1]>1) {FlowDirStart[1:(OutletPos[1]-1),(Transept+1):dimX] <- 7}
        FlowDirStart[OutletPos[1],(Transept+1):dimX] <- 1
        FlowDirStart[,1:Transept] <- 1
      } else if (FirstSide=="S"){
        Transept <- round(0.5*dimY)
        if (OutletPos[1]>1) {FlowDirStart[1:Transept,1:(OutletPos[1]-1)] <- 1}
        if (dimX>OutletPos[1]) {FlowDirStart[1:Transept,(OutletPos[1]+1):dimX] <- 5}
        FlowDirStart[1:Transept,OutletPos[1]] <- 3
        FlowDirStart[(Transept+1):dimY,] <- 3
      } else if (FirstSide=="W"){
        Transept <- round(0.5*dimX)
        if (dimY>OutletPos[1]) {FlowDirStart[(OutletPos[1]+1):dimY,1:Transept] <- 3}
        if (OutletPos[1]>1) {FlowDirStart[1:(OutletPos[1]-1),1:Transept] <- 7}
        FlowDirStart[OutletPos[1],1:Transept] <- 5
        FlowDirStart[,(Transept+1):dimX] <- 5
      }
      # impose drains perpendicular to outlets
      tmpX <- round(0.15*dimX); tmpY <- round(0.15*dimY)
      for (i in 1:N_outlet){
        if (OutletSide[i]=="N"){
          lowX <- max(1,OutletPos[i]-tmpX)
          uppX <- min(OutletPos[i]+tmpX,dimX)
          FlowDirStart[(dimY-tmpY+1):dimY,lowX:OutletPos[i]] <- 1
          FlowDirStart[(dimY-tmpY+1):dimY,OutletPos[i]:uppX] <- 5
        } else if (OutletSide[i]=="E") {
          lowY <- max(1,OutletPos[i]-tmpY)
          uppY <- min(OutletPos[i]+tmpY,dimY)
          FlowDirStart[lowY:OutletPos[i],(dimX-tmpX+1):dimX] <- 7
          FlowDirStart[OutletPos[i]:uppY,(dimX-tmpX+1):dimX] <- 3
        } else if (OutletSide[i]=="S") {
          lowX <- max(1,OutletPos[i]-tmpX)
          uppX <- min(OutletPos[i]+tmpX,dimX)
          FlowDirStart[1:tmpY,lowX:OutletPos[i]] <- 1
          FlowDirStart[1:tmpY,OutletPos[i]:uppX] <- 5
        } else if (OutletSide[i]=="W") {
          lowY <- max(1,OutletPos[i]-tmpY)
          uppY <- min(OutletPos[i]+tmpY,dimY)
          FlowDirStart[lowY:OutletPos[i],1:tmpX] <- 7
          FlowDirStart[OutletPos[i]:uppY,1:tmpX] <- 3
        }
      }
      for (i in 1:N_outlet){
        if (OutletSide[i]=="N"){
      FlowDirStart[(dimY-tmpY+1):dimY,OutletPos[i]] <- 7 
        } else if (OutletSide[i]=="E") {
      FlowDirStart[OutletPos[i],(dimX-tmpX+1):dimX] <- 1 
        } else if (OutletSide[i]=="S") {
      FlowDirStart[1:tmpY,OutletPos[i]] <- 3 
        } else if (OutletSide[i]=="W") {
      FlowDirStart[OutletPos[i],1:tmpX] <- 5 
        }
      }      
    } else if (type_initialstate=="V") {
      FirstSide <- OutletSide[1]
      if (FirstSide=="N"){
        for (i in 1:min(OutletPos[1],dimY)) {FlowDirStart[dimY+1-i,OutletPos[1]+1-i]=8}
        if (dimX>OutletPos[1]) {for (i in 1:min(dimX-OutletPos[1],dimY-1)) {FlowDirStart[dimY-i,OutletPos[1]+i]=6}}
        if (OutletPos[1]>1) {for (i in 1:min(OutletPos[1]-1,dimY)) {FlowDirStart[dimY+1-i,1:(OutletPos[1]-i)]=1}}
        if (dimX>OutletPos[1]) {for (i in 1:min(dimY,dimX-OutletPos[1])) {FlowDirStart[dimY+1-i,(OutletPos[1]+i):dimX]=5}}
        FlowDirStart[FlowDirStart==0]=7
      }
      if (FirstSide=="E"){
        for (i in 1:min(OutletPos[1],dimY)) {FlowDirStart[i,dimX-OutletPos[1]+i]=8}
        if (dimY>OutletPos[1]) {for (i in 1:min(dimY-OutletPos[1],dimX-1)) {FlowDirStart[OutletPos[1]+i,dimX-i]=2}}
        if (OutletPos[1]>1) {for (i in 1:min((OutletPos[1]-1),dimY)) {FlowDirStart[i,(dimX-OutletPos[1]+i+1):dimX]=7}}
        if (dimY>OutletPos[1]) {for (i in 1:min(dimY-OutletPos[1],dimX)) {FlowDirStart[(OutletPos[1]+i):dimY,dimX+1-i]=3}}
        FlowDirStart[FlowDirStart==0]=1
      }
      if (FirstSide=="S"){
        for (i in 1:min(OutletPos[1],dimY)) {FlowDirStart[i,OutletPos[1]+1-i]=2}
        if (dimX>OutletPos[1]) {for (i in 1:min(dimX-OutletPos[1],dimY-1)) {FlowDirStart[i+1,OutletPos[1]+i]=4}}
        if (OutletPos[1]>1) {for (i in 1:min(OutletPos[1]-1,dimY)) {FlowDirStart[i,1:(OutletPos[1]-i)]=1}}
        if (dimX>OutletPos[1]) {for (i in 1:min(dimX-OutletPos[1],dimY)) {FlowDirStart[i,(OutletPos[1]+i):dimX]=5}}
        FlowDirStart[FlowDirStart==0]=3
      }
      if (FirstSide=="W"){
        for (i in 1:min(OutletPos[1],dimY)) {FlowDirStart[i,OutletPos[1]-i+1]=6}
        if (dimY>OutletPos[1]) {for (i in 1:min(dimY-OutletPos[1],dimX-1)) {FlowDirStart[OutletPos[1]+i,i+1]=4}}
        if (OutletPos[1]>1) {for (i in 1:min(OutletPos[1]-1,dimY)) {FlowDirStart[i,1:(OutletPos[1]-i)]=7}}
        if (dimY>OutletPos[1]) {for (i in 1:min(dimY-OutletPos[1],dimX)) {FlowDirStart[(OutletPos[1]+i):dimY,i]=3}}
        FlowDirStart[FlowDirStart==0]=5
      }
      
      # impose drains with V shape
      # for (i in 1:N_outlet){
      #   if (OutletSide[i]=="N"){
      #     for (j in 1:min(OutletPos[i],dimY)) {FlowDirStart[dimY+1-j,OutletPos[i]+1-j]=8}
      #     if (dimX>OutletPos[i]) {for (j in 1:min(dimX-OutletPos[i],dimY-1)) {FlowDirStart[dimY-j,OutletPos[i]+j]=6}}
      #   } else if (OutletSide[i]=="E") {
      #     for (j in 1:OutletPos[i]) {FlowDirStart[j,dimX-OutletPos[i]+j]=8}
      #     if (dimY>OutletPos[i]) {for (j in 1:min(dimY-OutletPos[i],dimX-1)) {FlowDirStart[OutletPos[i]+j,dimX-j]=2}}
      #   } else if (OutletSide[i]=="S") {
      #     for (j in 1:min(OutletPos[i],dimY)) {FlowDirStart[j,OutletPos[i]+1-j]=2}
      #     if (dimX>OutletPos[i]) {for (j in 1:min(dimX-OutletPos[i],dimY-1)) {FlowDirStart[j+1,OutletPos[i]+j]=4}}
      #   } else if (OutletSide[i]=="W") {
      #     for (j in 1:OutletPos[i]) {FlowDirStart[j,OutletPos[i]-j+1]=6}
      #     if (dimY>OutletPos[i]) {for (j in 1:min(dimY-OutletPos[i],dimX-1)) {FlowDirStart[OutletPos[i]+j,j+1]=4}}
      #   }
      # }
      tmpX <- round(0.15*dimX); tmpY <- round(0.15*dimY)
      for (i in 1:N_outlet){
           if (OutletSide[i]=="N"){
             lowX <- max(1,OutletPos[i]-tmpX)
             uppX <- min(OutletPos[i]+tmpX,dimX)
             FlowDirStart[(dimY-tmpY):dimY,lowX:uppX] <- 0
             if (OutletPos[i]>1) {for (j in 1:(OutletPos[i]-lowX)) {FlowDirStart[dimY+1-j,lowX:(OutletPos[i]-j)]=1}}
             if (dimX>OutletPos[i]) {for (j in 1:(uppX-OutletPos[i])) {FlowDirStart[dimY+1-j,(OutletPos[i]+j):uppX]=5}}
             FlowDirStart[FlowDirStart==0]=7
             
           } else if (OutletSide[i]=="S"){
             lowX <- max(1,OutletPos[i]-tmpX)
             uppX <- min(OutletPos[i]+tmpX,dimX)
             FlowDirStart[1:(tmpY+1),lowX:uppX] <- 0
             if (OutletPos[i]>1) {for (j in 1:(OutletPos[i]-lowX)) {FlowDirStart[j,lowX:(OutletPos[i]-j)]=1}}
             if (dimX>OutletPos[i]) {for (j in 1:(uppX-OutletPos[i])) {FlowDirStart[j,(OutletPos[i]+j):uppX]=5}}
             FlowDirStart[FlowDirStart==0]=3
             
           } else if (OutletSide[i]=="E"){
             lowY <- max(1,OutletPos[i]-tmpY)
             uppY <- min(OutletPos[i]+tmpY,dimY)
             FlowDirStart[lowY:uppY,(dimX-tmpX):dimX] <- 0
             if (OutletPos[i]>1) {for (j in 1:(OutletPos[i]-lowY)) {FlowDirStart[lowY-1+j,(dimX-OutletPos[i]+lowY+j-1):dimX]=7}}
             if (dimY>OutletPos[i]) {for (j in 1:(uppY-OutletPos[i])) {FlowDirStart[OutletPos[i]+j,(dimX+1-j):dimX]=3}}
             FlowDirStart[FlowDirStart==0]=1
           }  else if (OutletSide[i]=="W"){
             lowY <- max(1,OutletPos[i]-tmpY)
             uppY <- min(OutletPos[i]+tmpY,dimY)
             FlowDirStart[lowY:uppY,1:(tmpX+1)] <- 0
             if (OutletPos[i]>1) {for (j in 1:(OutletPos[i]-lowY)) {FlowDirStart[lowY-1+j,1:(OutletPos[i]-lowY-j+1)]=7}}
             if (dimY>OutletPos[i]) {for (j in 1:(uppY-OutletPos[i])) {FlowDirStart[OutletPos[i]+j,1:j]=3}}
             FlowDirStart[FlowDirStart==0]=5
           }
        }
      for (i in 1:N_outlet){
        if (OutletSide[i]=="N"){
          lowX <- max(1,OutletPos[i]-tmpX)
          uppX <- min(OutletPos[i]+tmpX,dimX)
          if (OutletPos[i]>1) {for (j in 1:(OutletPos[i]-lowX+1)) {FlowDirStart[dimY+1-j,OutletPos[i]+1-j]=8}}
          if (dimX>OutletPos[i]) {for (j in 1:(uppX-OutletPos[i])) {FlowDirStart[dimY-j,OutletPos[i]+j]=6}}
          
        } else if (OutletSide[i]=="S"){
          lowX <- max(1,OutletPos[i]-tmpX)
          uppX <- min(OutletPos[i]+tmpX,dimX)
          if (OutletPos[i]>1) {for (j in 1:(OutletPos[i]-lowX+1)) {FlowDirStart[j,OutletPos[i]+1-j]=2}}
          if (dimX>OutletPos[i]) {for (j in 1:(uppX-OutletPos[i])) {FlowDirStart[j+1,OutletPos[i]+j]=4}}
          
        } else if (OutletSide[i]=="E"){
          lowY <- max(1,OutletPos[i]-tmpY)
          uppY <- min(OutletPos[i]+tmpY,dimY)
          if (OutletPos[i]>1) {for (j in 1:(OutletPos[i]-lowY)) {FlowDirStart[lowY-1+j,dimX-OutletPos[i]+lowY+j-1]=8}}
          if (dimY>OutletPos[i]) {for (j in 1:(uppY-OutletPos[i])) {FlowDirStart[OutletPos[i]+j,dimX-j]=2}}
        }  else if (OutletSide[i]=="W"){
          lowY <- max(1,OutletPos[i]-tmpY)
          uppY <- min(OutletPos[i]+tmpY,dimY)
          if (OutletPos[i]>1) {for (j in 1:(OutletPos[i]-lowY)) {FlowDirStart[lowY-1+j,OutletPos[i]-lowY-j+2]=6}}
          if (dimY>OutletPos[i]) {for (j in 1:(uppY-OutletPos[i])) {FlowDirStart[OutletPos[i]+j,j+1]=4}}
        }
      }
      
    } else {
      stop("Invalid initial state")
    }
  }
  return(FlowDirStart)
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
  
  return(OutList)
}


