
draw_elev2D_OCN <- function(OCN,
                              colPalette=terrain.colors(1000,alpha=1),
                              addLegend=TRUE){
  
  if (!("Z" %in% names(OCN$FD))){
    stop('Missing fields in OCN. You should run landscape_OCN prior to draw_elev2D_OCN.')
  }

  # plot elevation map
  #old.par <- par(no.readonly =TRUE) 
  #on.exit(par(old.par))
  #par(bty="n")
  
  if (is.null(OCN$xllcorner)){xllcorner <- min(OCN$FD$X)[1]} else {xllcorner <- OCN$xllcorner}
  if (is.null(OCN$yllcorner)){yllcorner <- min(OCN$FD$Y)[1]} else {yllcorner <- OCN$yllcorner}
  
  if (OCN$FD$nNodes < OCN$dimX*OCN$dimY){
    Zmat <- matrix(NaN,OCN$dimY,OCN$dimX)
    Zmat[OCN$FD$toDEM] <- OCN$FD$Z
    Zmat <- Zmat[seq(OCN$dimY,1,-1),]
  } else {Zmat <- matrix(data=OCN$FD$Z,nrow=OCN$dimY,ncol=OCN$dimX)}
  
  if (addLegend==TRUE){
  image.plot(seq(xllcorner,xllcorner+(OCN$dimX-1)*OCN$cellsize,OCN$cellsize),
             seq(yllcorner,yllcorner+(OCN$dimY-1)*OCN$cellsize,OCN$cellsize),
             t(Zmat),col=colPalette,xlab=" ",ylab=" ",axes=FALSE,asp=1)
  } else {
    image(seq(xllcorner,xllcorner+(OCN$dimX-1)*OCN$cellsize,OCN$cellsize),
          seq(yllcorner,yllcorner+(OCN$dimY-1)*OCN$cellsize,OCN$cellsize),
               t(Zmat),col=colPalette,xlab=" ",ylab=" ",axes=FALSE,asp=1)
  }
  invisible()
}