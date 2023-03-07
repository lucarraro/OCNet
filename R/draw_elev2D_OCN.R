
draw_elev2D_OCN <- function(OCN,
                            colPalette=terrain.colors(1000,alpha=1),
                            addLegend=TRUE,
                            drawRiver=FALSE,
                            thrADraw=0.002*OCN$FD$nNodes*OCN$cellsize^2,
                            riverColor="#00BFFF",
                            min_lwd=0.5,
                            max_lwd=5,
                            args_imagePlot=list()){
  
  if (!("Z" %in% names(OCN$FD))){
    stop('Missing fields in OCN. You should run landscape_OCN prior to draw_elev2D_OCN.')
  }
  
  # plot elevation map
  #old.par <- par(no.readonly =TRUE) 
  #on.exit(par(old.par))
  #par(bty="n")
  
  if (is.null(OCN$xllcorner)){xllcorner <- min(OCN$FD$X)[1]} else {xllcorner <- OCN$xllcorner}
  if (is.null(OCN$yllcorner)){yllcorner <- min(OCN$FD$Y)[1]} else {yllcorner <- OCN$yllcorner}
  
  Xvec <- seq(xllcorner,xllcorner+(OCN$dimX-1)*OCN$cellsize,OCN$cellsize)
  Yvec <- seq(yllcorner,yllcorner+(OCN$dimY-1)*OCN$cellsize,OCN$cellsize)
  
  if (OCN$FD$nNodes < OCN$dimX*OCN$dimY){
    if (isTRUE(OCN$typeInitialState=="custom")){
      Zmat <- matrix(NaN,OCN$dimY,OCN$dimX)
      Zmat[OCN$FD$toDEM] <- OCN$FD$Z
      Zmat <- Zmat[seq(OCN$dimY,1,-1), ]
    } else { # real river
      Zmat <- matrix(NaN,OCN$dimX,OCN$dimY)
      Zmat[OCN$FD$toDEM] <- OCN$FD$Z
      Zmat <- Zmat[,seq(OCN$dimY,1,-1)]
      Zmat <- t(Zmat)
    }
  } else {Zmat <- matrix(data=OCN$FD$Z,nrow=OCN$dimY,ncol=OCN$dimX)}
  
  if(is.null(args_imagePlot$col)){args_imagePlot$col <- colPalette}
  if(is.null(args_imagePlot$zlim)){args_imagePlot$zlim <- range(Zmat,na.rm=T)}
  args_imagePlot$x <- seq(xllcorner,xllcorner+(OCN$dimX-1)*OCN$cellsize,OCN$cellsize)
  args_imagePlot$y <- seq(yllcorner,yllcorner+(OCN$dimY-1)*OCN$cellsize,OCN$cellsize)
  args_imagePlot$z <- t(Zmat)
  if(is.null(args_imagePlot$xlab)){args_imagePlot$xlab <- " "} 
  if(is.null(args_imagePlot$ylab)){args_imagePlot$ylab <- " "}
  if(is.null(args_imagePlot$axes)){args_imagePlot$axes <- FALSE} 
  if(is.null(args_imagePlot$asp)){args_imagePlot$asp <- 1}
  if(is.null(args_imagePlot$xlim)){args_imagePlot$xlim=range(OCN$FD$X)}
  if(is.null(args_imagePlot$ylim)){args_imagePlot$ylim=range(OCN$FD$Y)}
  
  if (addLegend==TRUE){
    do.call(imagePlot, args_imagePlot)
  } else {do.call(image, args_imagePlot)}
  
  if (drawRiver==TRUE){
    AvailableNodes <- setdiff(1:OCN$FD$nNodes,OCN$FD$outlet)
    for (i in AvailableNodes){
      if (OCN$FD$A[i]>=thrADraw  & abs(OCN$FD$X[i]-OCN$FD$X[OCN$FD$downNode[i]]) <= 1.001*OCN$cellsize & abs(OCN$FD$Y[i]-OCN$FD$Y[OCN$FD$downNode[i]]) <= 1.001*OCN$cellsize) {
        lines(c(OCN$FD$X[i],OCN$FD$X[OCN$FD$downNode[i]]),c(OCN$FD$Y[i],OCN$FD$Y[OCN$FD$downNode[i]]),
              lwd=min_lwd+(max_lwd-min_lwd)*(OCN$FD$A[i]/(OCN$FD$nNodes*OCN$cellsize^2))^0.5,col=riverColor)}}
  }
  
  invisible()
}