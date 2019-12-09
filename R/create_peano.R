
create_peano <- function(nIterPeano,
                         outletPos="NE",
                         xllcorner=1,
                         yllcorner=1,
                         cellsize=1){
  
  if (!(outletPos %in% c("NE","NW","SE","SW"))){
    stop('Invalid outletPos')}
  
  # create Peano flow direction matrix
  FlowDir_peano <- matrix(c(8,7,1,8),nrow=2,ncol=2,byrow=TRUE)
  
  for (i in 1:nIterPeano){
    FlowDir_peano <- rbind(cbind(FlowDir_peano,flip_lr_peano(FlowDir_peano)),
                           cbind(flip_ud_peano(FlowDir_peano),FlowDir_peano))
  }
  dim <- 2^(nIterPeano+1)
  FlowDir_peano[dim,dim] <- 0
  Nnodes <- dim^2
  Outlet <- Nnodes
  
  # flip flow direction matrix if outlet position is not NE
  if (outletPos=="SE"){
    FlowDir_peano <- flip_ud_peano(FlowDir_peano)
    FlowDir_peano[1,dim] <- 0
    Outlet <- Nnodes-dim+1}
  if (outletPos=="NW"){
    FlowDir_peano <- flip_lr_peano(FlowDir_peano)
    FlowDir_peano[dim,1] <- 0
    Outlet <- dim}
  if (outletPos=="SW"){
    FlowDir_peano <- flip_ud_peano(flip_lr_peano(FlowDir_peano))
    FlowDir_peano[1,1] <- 0
    FlowDir_peano[1,dim] <- 6
    Outlet <- 1}
  
  # find list of possible neighbouring pixels
  
  movement <- matrix(c(0,-1,-1,-1,0,1,1,1,1,1,0,-1,-1,-1,0,1),nrow=2,byrow=TRUE)
  
  # Adjacency matrix for the initial state
  cont_node <- 0
  W <- spam(0,Nnodes,Nnodes)
  ind <- matrix(0,Nnodes,2)
  DownNode <- numeric(Nnodes)
  for (cc in 1:dim) {
    for (rr in 1:dim) {
      cont_node <- cont_node + 1
      dir <- FlowDir_peano[rr,cc]
      if (dir>0) {
        DownPixel <- c(rr,cc)+movement[,dir]
        #AD[cont_node,(DownPixel[2]-1)*dim+DownPixel[1]] <- 1
        DownNode[cont_node] <- (DownPixel[2]-1)*dim+DownPixel[1]
        ind[cont_node,] <- c(cont_node,(DownPixel[2]-1)*dim+DownPixel[1])
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
  
  pl <- initial_permutation(DownNode)
  
  pas <- permuteAddSolve(Wt,pl$perm,1L, numeric(Nnodes),0.5)#*(cellsize^2)
  A <- pas[[1]]
  A <- A[invPerm(pl$perm)]
  
  # define identity matrix and vector of ones (needed to calculate A)
  #Imat <- sparseMatrix(i=1:Nnodes,j=1:Nnodes)
  #Ones <- rep(1,Nnodes)
  # calculate vector of contributing area
  #A <- solve(Imat-t(AD),Ones)
  
  # generate X, Y, NeighbouringNodes (for compatibility with function create_ocn)
  X <- rep((xllcorner*cellsize):cellsize:(dim*cellsize), each = dim)  
  Y <- rep((yllcorner*cellsize):cellsize:(dim*cellsize),dim)
  
  
  # find list of possible neighbouring pixels
  movement <- matrix(c(0,-1,-1,-1,0,1,1,1,1,1,0,-1,-1,-1,0,1),nrow=2,byrow=TRUE)
  NeighbouringNodes <- vector("list", Nnodes)
  cont_node <- 0
  for (cc in 1:dim) {
    for (rr in 1:dim) {
      cont_node <- cont_node + 1
      neigh_r <- rep(rr,8)+movement[1,]
      neigh_c <- rep(cc,8)+movement[2,]
      NotAboundary <- neigh_r>0 & neigh_r<=dim & neigh_c>0 & neigh_c<=dim 
      NeighbouringNodes[[cont_node]] <- neigh_r[NotAboundary] + (neigh_c[NotAboundary]-1)*dim
    }
  } 
  
  # write results 
  ## create list for export
  FD <- list(A=A,W=W,downNode=DownNode,X=X,Y=Y,nNodes=Nnodes,outlet=Outlet)
  peano <- list(FD=FD,dimX=dim,dimY=dim,cellsize=cellsize,nOutlet=1,periodicBoundaries=FALSE,expEnergy=0.5)
  # note: for peano networks, ExpEnergy only has the meaning of quantity related to the exponent of the slope-area relationship
  # s \propto A^(ExpEnergy-1)
  
  invisible(peano)   
}

# auxiliary functions
flip_lr_peano <- function(mat){
  # flip matrix 
  dim <- ncol(mat)
  mat <- mat[,dim:1]
  # mirror
  mat[mat==1] <- 50
  mat[mat==2] <- 40
  mat[mat==3] <- 30
  mat[mat==4] <- 20
  mat[mat==5] <- 10
  mat[mat==6] <- 80
  mat[mat==7] <- 70
  mat[mat==8] <- 60
  mat <- mat/10 
  mat[dim,1] <- 7
  invisible(mat)
}

flip_ud_peano <- function(mat){
  dim <- ncol(mat)
  mat <- mat[dim:1,]
  # mirror
  mat[mat==1] <- 10
  mat[mat==2] <- 80
  mat[mat==3] <- 70
  mat[mat==4] <- 60
  mat[mat==5] <- 50
  mat[mat==6] <- 40
  mat[mat==7] <- 30
  mat[mat==8] <- 20
  mat <- mat/10 
  mat[1,dim] <- 1
  invisible(mat)
}



