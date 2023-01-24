
draw_subcatchments_OCN <- function(OCN,
                                   drawRiver = TRUE,
                                   colPalette = NULL){
  
  if (!("SC" %in% names(OCN))){
    stop('Missing fields in OCN. You should run aggregate_OCN prior to draw_subcatchments_OCN.')
  }
  
  
  ## Greedy algorithm for coloring subcatchment map
  ## Create list of nodes for the greedy algorithm
  W <- OCN$SC$W
  ListNodes <- numeric(OCN$SC$nNodes)
  for (i in 1:(OCN$SC$nNodes)){
    degree <- colSums(W)
    degree <- degree[degree>0]
    if (length(degree)>0){
      MinDegree <- which(colSums(W)==min(degree))
      NodeRemoved <- MinDegree[1]
      ListNodes[OCN$SC$nNodes+1-i] <- NodeRemoved
      W[,NodeRemoved] <- 0
      W[NodeRemoved,] <- 0
    } else {
      MissingNode <- setdiff(1:OCN$SC$nNodes,ListNodes)[1]
      ListNodes[OCN$SC$nNodes+1-i] <- MissingNode
    }
  }
  
  
  ## Attribute colors
  ColorList <- 1
  ColorID <- numeric(OCN$SC$nNodes)
  for (i in 1:OCN$SC$nNodes){
    connected_nodes <- which(OCN$SC$W[ListNodes[i],1:OCN$SC$nNodes]==1)
    k <- setdiff(ColorList,ColorID[connected_nodes])[1]
    if (is.na(k)){
      ColorList <- c(ColorList,max(ColorList)+1)
      k <- max(ColorList)
    }
    ColorID[ListNodes[i]] <- k
  }
  
  ## plot subcatchment map
  kol <- numeric(OCN$FD$nNodes)
  for (k in 1:OCN$SC$nNodes){
    kol[OCN$SC$toFD[[k]]] <- ColorID[k]
  }
  
  if (OCN$FD$nNodes < OCN$dimX*OCN$dimY){
    Color_SC <- matrix(NaN,OCN$dimY,OCN$dimX)
    Color_SC[OCN$FD$toDEM] <- kol
    Color_SC <- Color_SC[seq(OCN$dimY,1,-1),]
  } else {
    Color_SC <- matrix(data=kol,nrow=OCN$dimY,ncol=OCN$dimX)
  }
  
  
  if (is.null(colPalette)){
  colPalette <- c("#009900", # green
                  "#FFFF00", # yellow
                  "#FF9900", # orange
                  "#FF0000", # red
                  "#FF00FF", # fuchsia
                  "#9900CC", # violet
                  "#555555", # grey 30%
                  "#BBBBBB") # grey 70%
  
  colPalette <- colPalette[ColorList]
  } else if (typeof(colPalette)=="closure") {
    colPalette <- colPalette(length(ColorList))
  } else if (typeof(colPalette)=="character") {
    colPalette <- colPalette[1:length(ColorList)]
  }
  
  if (is.null(OCN$xllcorner)){xllcorner <- min(OCN$FD$X)[1]} else {xllcorner <- OCN$xllcorner}
  if (is.null(OCN$yllcorner)){yllcorner <- min(OCN$FD$Y)[1]} else {yllcorner <- OCN$yllcorner}
  
  Xvec <- seq(xllcorner,xllcorner+(OCN$dimX-1)*OCN$cellsize,OCN$cellsize)
  Yvec <- seq(yllcorner,yllcorner+(OCN$dimY-1)*OCN$cellsize,OCN$cellsize)
  
  #old.par <- par(no.readonly = TRUE)
  #on.exit(par(old.par))
  #par(bty="n")
  image(Xvec, Yvec,
        t(Color_SC),col=colPalette,xlab=" ",ylab=" ",asp=1, axes=FALSE) #
  # attributing colors in reverse order should increase overall contrast
  
  if (drawRiver==TRUE){
    ## plot OCN
    AvailableNodes <- setdiff(1:OCN$FD$nNodes,OCN$FD$outlet)
    #points(OCN$FD$X[OCN$FD$outlet],OCN$FD$Y[OCN$FD$outlet],pch=22,col="#000000",bg="#000000")

    for (i in AvailableNodes){
      if (OCN$FD$A[i]>=OCN$thrA  & abs(OCN$FD$X[i]-OCN$FD$X[OCN$FD$downNode[i]])<=OCN$cellsize & abs(OCN$FD$Y[i]-OCN$FD$Y[OCN$FD$downNode[i]])<=OCN$cellsize  ) {
        lines(c(OCN$FD$X[i],OCN$FD$X[OCN$FD$downNode[i]]),c(OCN$FD$Y[i],OCN$FD$Y[OCN$FD$downNode[i]]),lwd=0.5+4.5*(OCN$FD$A[i]/(OCN$FD$nNodes*OCN$cellsize^2))^0.5,col="black")}
    }
  }
  invisible()
}