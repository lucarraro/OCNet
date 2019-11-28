
draw_subcatchments_OCN <- function(OCN,
                                   draw_river = TRUE,
                                   ColPalette = NULL){
  
  if (!("SC" %in% names(OCN))){
    stop('Missing fields in OCN. You should run aggregate_OCN prior to draw_subcatchments_OCN.')
  }
  
  ## Greedy algorithm for coloring subcatchment map
  ## Create list of nodes for the greedy algorithm
  W <- OCN$SC$W
  ListNodes <- numeric(OCN$SC$Nnodes)
  for (i in 1:(OCN$SC$Nnodes)){
    degree <- colSums(W)
    degree <- degree[degree>0]
    if (length(degree)>0){
      MinDegree <- which(colSums(W)==min(degree))
      NodeRemoved <- MinDegree[1]
      ListNodes[OCN$SC$Nnodes+1-i] <- NodeRemoved
      W[,NodeRemoved] <- 0
      W[NodeRemoved,] <- 0
    } else {
      MissingNode <- setdiff(1:OCN$SC$Nnodes,ListNodes)[1]
      ListNodes[OCN$SC$Nnodes+1-i] <- MissingNode
    }
  }
  
  
  ## Attribute colors
  ColorList <- 1
  ColorID <- numeric(OCN$SC$Nnodes)
  for (i in 1:OCN$SC$Nnodes){
    connected_nodes <- which(OCN$SC$W[ListNodes[i],1:OCN$SC$Nnodes]==1)
    k <- setdiff(ColorList,ColorID[connected_nodes])[1]
    if (is.na(k)){
      ColorList <- c(ColorList,max(ColorList)+1)
      k <- max(ColorList)
    }
    ColorID[ListNodes[i]] <- k
  }
  
  ## plot subcatchment map
  Color_SC <- matrix(data=OCN$FD$to_SC,nrow=OCN$dimY,ncol=OCN$dimX)
  for (k in 1:OCN$SC$Nnodes){
    Color_SC[OCN$SC$to_FD[[k]]] <- ColorID[k]
  }
  
  
  if (is.null(ColPalette)){
  ColPalette <- c("#009900", # green
                  "#FFFF00", # yellow
                  "#FF9900", # orange
                  "#FF0000", # red
                  "#FF00FF", # fuchsia
                  "#9900CC", # violet
                  "#555555", # grey 30%
                  "#BBBBBB") # grey 70%
  
  ColPalette <- ColPalette[ColorList]
  } else if (typeof(ColPalette)=="closure") {
    ColPalette <- ColPalette(length(ColorList))
  } else if (typeof(ColPalette)=="character") {
    ColPalette <- ColPalette[1:length(ColorList)]
  }
  
  par(pty="s",bty="n",mar=c(1,1,1,1))
  image(seq(min(OCN$FD$X),max(OCN$FD$X),OCN$cellsize),
        seq(min(OCN$FD$Y),max(OCN$FD$Y),OCN$cellsize),
        t(Color_SC),col=ColPalette,xlab=" ",ylab=" ",asp=1,axes=FALSE)
  # attributing colors in reverse order should increase overall contrast
  
  if (draw_river==TRUE){
  ## plot OCN
  AvailableNodes <- setdiff(1:OCN$FD$Nnodes,OCN$FD$Outlet)
  #points(OCN$FD$X[OCN$FD$Outlet],OCN$FD$Y[OCN$FD$Outlet],pch=22,col="#000000",bg="#000000")

  for (i in AvailableNodes){
    if (OCN$FD$A[i]>=OCN$A_thr  & abs(OCN$FD$X[i]-OCN$FD$X[OCN$FD$DownNode[i]])<=OCN$cellsize & abs(OCN$FD$Y[i]-OCN$FD$Y[OCN$FD$DownNode[i]])<=OCN$cellsize  ) {
      lines(c(OCN$FD$X[i],OCN$FD$X[OCN$FD$DownNode[i]]),c(OCN$FD$Y[i],OCN$FD$Y[OCN$FD$DownNode[i]]),lwd=0.5+4.5*(OCN$FD$A[i]/(OCN$FD$Nnodes*OCN$cellsize^2))^0.5,col="black")}
  }
  }
}