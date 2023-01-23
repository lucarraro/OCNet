create_general_contour_OCN <- function(flowDirStart,
                       expEnergy=0.5, 
                       cellsize=1,
                       xllcorner=0.5*cellsize,
                       yllcorner=0.5*cellsize,
                       nIter=NULL,
                       nUpdates=50,
                       initialNoCoolingPhase=0,
                       coolingRate=1,
                       showIntermediatePlots=FALSE,
                       thrADraw=NULL,
                       easyDraw=NULL,
                       saveEnergy=FALSE,
                       saveExitFlag=FALSE,
                       displayUpdates=1){
  
  periodicBoundaries <- FALSE
  typeInitialState <- "custom"
  
  Nnodes_FD <- sum(!is.na(flowDirStart))
  if (is.null(thrADraw)){thrADraw <- 0.002*Nnodes_FD*cellsize^2}
  if (is.null(nIter)){nIter <- 40*Nnodes_FD}
  
  if (Nnodes_FD > 1000) {
    dim <- sqrt(Nnodes_FD)
    estTime <- (7.988e-1 - 1.266e-1*dim + 7.198e-3*dim^2 - 6.807e-5*dim^3 +  1.372e-6*dim^4) / (30*Nnodes_FD) * nIter
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
      message("Note that the above estimate is only based on the choice of parameter nIter, and not on processor performance.\n", appendLF = FALSE)
      message("\n", appendLF = FALSE)
    }
  }
  
  t0 <- Sys.time()
  if (displayUpdates==2){message('Initializing...\n', appendLF = FALSE)}


  if (is.null(easyDraw)){
    if (Nnodes_FD>4e4) {
      easyDraw=TRUE
    } else {easyDraw=FALSE}
  }
  
  if (sum(flowDirStart==0,na.rm=T)==0){
    stop("flowDirStart has no imposed outlet(s).")
  }

t0 <- Sys.time()

Nnodes_FD <- sum(!is.na(flowDirStart))
dimX <- ncol(flowDirStart)
dimY <- nrow(flowDirStart)
Nnodes_DEM = dimX*dimY
FD_to_DEM <- which(!is.nan(flowDirStart))
DEM_to_FD <- numeric(Nnodes_DEM)
DEM_to_FD[which(!is.nan(flowDirStart))] <- 1:Nnodes_FD

Xvec <- seq(xllcorner, xllcorner+(dimX-1)*cellsize, cellsize)
Yvec <- seq(xllcorner, yllcorner+(dimY-1)*cellsize, cellsize)
X <- Y <- numeric(Nnodes_FD)
for (i in 1:Nnodes_FD){
  indDEM <- FD_to_DEM[i]
  rr <- indDEM %% dimY + ((indDEM %% dimY)==0)*dimY
  cc <- (indDEM - rr)/dimY + 1
  X[i] <- Xvec[cc]
  Y[i] <- Yvec[dimY-rr]
}

# find list of possible neighbouring pixels
movement <- matrix(c(0,-1,-1,-1,0,1,1,1,1,1,0,-1,-1,-1,0,1),nrow=2,byrow=TRUE)
NeighbouringNodes <- vector("list", Nnodes_DEM)
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

NeighbouringNodes_FD <- vector("list", Nnodes_FD)
for (i in 1:Nnodes_FD){
  indDEM <- FD_to_DEM[i]
  tmp <- DEM_to_FD[NeighbouringNodes[[indDEM]]]
  NeighbouringNodes_FD[[i]] <- tmp[tmp != 0]
}

# DEM level: pixels of the extracted DEM
# FD level: every node is a catchment node

W <- spam(0,Nnodes_FD,Nnodes_FD)
ind <- matrix(0,Nnodes_FD,2)
DownNode <- numeric(Nnodes_FD) 

for (i in 1:Nnodes_FD){
  indDEM <- FD_to_DEM[i]
  dir <- flowDirStart[indDEM]
  if (dir > 0){
    rr <- indDEM %% dimY + ((indDEM %% dimY)==0)*dimY
    cc <- (indDEM - rr)/dimY + 1
    DownPixel <- c(rr,cc)+movement[,dir]
    ind[i,] <- c(i, DEM_to_FD[(DownPixel[2]-1)*dimY+DownPixel[1]])
    DownNode[i] <- DEM_to_FD[(DownPixel[2]-1)*dimY+DownPixel[1]]
  }
}
ind <- ind[-which(ind[,1]==0),]
W[ind] <- 1

OutletPixel <- which(DownNode==0)
nOutlet <- length(OutletPixel)

# patch to correct for initial 0
newW <- new("spam")
slot(newW, "entries", check = FALSE) <- W@entries[-1]
slot(newW, "colindices", check = FALSE) <- W@colindices[-1]
slot(newW, "rowpointers", check = FALSE) <- c(W@rowpointers[1],W@rowpointers[-1]-W@rowpointers[1])
slot(newW, "dimension", check = FALSE) <- W@dimension 

Wt <- t(newW)
rm(W,newW)

pl <- initial_permutation(DownNode)  # calculate permutation vector
pas <- permuteAddSolve(Wt, pl$perm, 1L, numeric(Nnodes_FD), 0.5)
A <- pas[[1]]
A <- A[invPerm(pl$perm)]

Amask <- flowDirStart
Amask[FD_to_DEM] <- A

# initialize energy and temperatures
Energy <- numeric(nIter)
Energy_0 <- pas[[2]]
Energy[1] <- Energy_0
Temperature <- c(Energy[1]+numeric(initialNoCoolingPhase*nIter),Energy[1]*exp(-coolingRate*(1:(nIter-initialNoCoolingPhase*nIter))/(Nnodes_FD)))
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
  catch <- numeric(Nnodes_FD)
  for (o in 1:nOutlet){
    catch[(resort[OutletPixel[o]]-AA[OutletPixel[o]]/cellsize^2+1):resort[OutletPixel[o]]] <- o
  }
  old.par <- par(bty="n")
  on.exit(par(old.par))
  plot(c(min(X),max(X)),c(min(Y),max(Y)),main=sprintf('OCN (initial state)'),
       type="n",asp=1,axes=FALSE,xlab="",ylab="") # 
  points(X[OutletPixel],Y[OutletPixel],pch=15,col=rnbw[catch[resort[OutletPixel]]])
  if (easyDraw == FALSE){
    for (i in AvailableNodes){
      if (AA[i]<=thrADraw & abs(X[i]-X[DownNode[i]])<=1.01*cellsize & abs(Y[i]-Y[DownNode[i]])<=1.01*cellsize ) {
        lines(c(X[i],X[DownNode[i]]),c(Y[i],Y[DownNode[i]]),lwd=0.5,col="#E0E0E0")}
    }
  }
  for (i in 1:Nnodes_FD){
    if (!(i %in% OutletPixel)){
      if (AA[i]>thrADraw & abs(X[i]-X[DownNode[i]])<=1.01*cellsize & abs(Y[i]-Y[DownNode[i]])<=1.01*cellsize ) {
        lines(c(X[i],X[DownNode[i]]),c(Y[i],Y[DownNode[i]]),lwd=0.5+4.5*(AA[i]/Nnodes_FD/cellsize^2)^0.5,col=rnbw[catch[resort[i]]])}
    }}
}

if (displayUpdates == 2){
  message(sprintf('Initialization completed. Elapsed time is %.2f s \n',difftime(Sys.time(),t0,units='secs')), appendLF = FALSE) 
  message('Search algorithm has started.\n', appendLF = FALSE)
  message('\n', appendLF = FALSE)}

pl <- as.integer(pl$perm)
flag <- 0

AvailableNodes <- setdiff(1:Nnodes_FD, OutletPixel)

# simulated annealing algorithm
if (nIter > 1){
  for (iter in 2:nIter) {
    t2 <- Sys.time() 
    # pick random node (excluding the outlet)
    Energy_new <- 100*Energy_0
    node <- sample(AvailableNodes,1)
    
    # change downstream connection from chosen node
    tmp <- NeighbouringNodes_FD[[node]][NeighbouringNodes_FD[[node]]!=DownNode[node]]
    if (length(tmp)==1){tmp <- c(tmp,tmp)} # patch for sample's surprise
    if (length(tmp)==0){tmp <- c(node,node)} # patch for no alternative directions
    down_new <- sample(tmp,1) # sample one node from list of neighbouring nodes, excluding the one that was previously connected to node
    
    ind1 <- (X[NeighbouringNodes_FD[[node]]]==X[node] & Y[NeighbouringNodes_FD[[node]]]==Y[down_new])
    Node1 <- NeighbouringNodes_FD[[node]][ind1]
    
    ind2 <- (X[NeighbouringNodes_FD[[node]]]==X[down_new] & Y[NeighbouringNodes_FD[[node]]]==Y[node])
    Node2 <- NeighbouringNodes_FD[[node]][ind2]
    
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
        catch <- numeric(Nnodes_FD)
        for (o in 1:nOutlet){
          catch[(resort[OutletPixel[o]]-AA[OutletPixel[o]]/cellsize^2+1):resort[OutletPixel[o]]] <- o
        }
        plot(c(min(X),max(X)),c(min(Y),max(Y)),type="n",main=sprintf('OCN (%.1f%% completed)',iter/nIter*100),
             asp=1,axes=FALSE,xlab=" ",ylab=" ") #
        points(X[OutletPixel],Y[OutletPixel],pch=15,col=rnbw[catch[resort[OutletPixel]]])
        if (easyDraw == FALSE){
          for (i in AvailableNodes){
            if (AA[i]<=(thrADraw)  & abs(X[i]-X[DownNode[i]])<=1.01*cellsize & abs(Y[i]-Y[DownNode[i]])<=1.01*cellsize ) {
              lines(c(X[i],X[DownNode[i]]),c(Y[i],Y[DownNode[i]]),lwd=0.5,col="#E0E0E0")}
          }
        }
        for (i in 1:Nnodes_FD){
          if (!(i %in% OutletPixel)){
            if (AA[i]>(thrADraw)  & abs(X[i]-X[DownNode[i]])<=1.01*cellsize & abs(Y[i]-Y[DownNode[i]])<=1.01*cellsize ) {
              lines(c(X[i],X[DownNode[i]]),c(Y[i],Y[DownNode[i]]),lwd=0.5+4.5*(AA[i]/Nnodes_FD/cellsize^2)^0.5,col=rnbw[catch[resort[i]]])}
          }}
      }}
  }
  t3<-Sys.time()
  savetime[iter] <- t3-t2
  ExitFlag[iter] <- flag
  
  if (sum(A[OutletPixel])  != Nnodes_FD){
    stop('Error: sum(A[OutletPixel]) is not equal to the total lattice area')
  }
}

W <- t(Wt)

FD <- list(A=A*cellsize^2,W=W,downNode=DownNode,X=X,Y=Y,nNodes=Nnodes_FD,outlet=OutletPixel,perm=pl,toDEM=FD_to_DEM)
OCN <- list(FD=FD,dimX=dimX,dimY=dimY,cellsize=cellsize,nOutlet=nOutlet,periodicBoundaries=periodicBoundaries,
            expEnergy=expEnergy,coolingRate=coolingRate,typeInitialState=typeInitialState,nIter=nIter,initialNoCoolingPhase=initialNoCoolingPhase,
            energyInit=Energy_0)

if (saveEnergy==TRUE) {OCN[["energy"]] <- Energy}
if (saveExitFlag==TRUE) {OCN[["exitFlag"]] <- ExitFlag}

invisible(OCN)
}