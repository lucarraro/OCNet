
draw_elev2D_OCN <- function(OCN,
                              colPalette=terrain.colors(1000,alpha=1),
                              addLegend=TRUE){
  
  if (!("Z" %in% names(OCN$FD))){
    stop('Missing fields in OCN. You should run landscape_OCN prior to draw_elev2D_OCN.')
  }
  if (isTRUE(OCN$typeInitialState=="custom")){
    stop('draw_elev2D_OCN is not currently implemented for OCNs created via create_general_contour_OCN')
  }
  
  # plot elevation map
  #old.par <- par(no.readonly =TRUE) 
  #on.exit(par(old.par))
  #par(bty="n")
  
  if (is.null(OCN$xllcorner)){xllcorner <- min(OCN$FD$x)} else {xllcorner <- OCN$xllcorner}
  if (is.null(OCN$yllcorner)){yllcorner <- min(OCN$FD$Y)} else {yllcorner <- OCN$yllcorner}
  
  Zmat<-matrix(data=OCN$FD$Z,nrow=OCN$dimY,ncol=OCN$dimX)
  if (addLegend==TRUE){
  image.plot(seq(xllcorner,xllcorner+OCN$dimX*OCN$cellsize,OCN$cellsize),
             seq(yllcorner,yllcorner+OCN$dimY*OCN$cellsize,OCN$cellsize),
             t(Zmat),col=colPalette,xlab=" ",ylab=" ",axes=FALSE,asp=1)
  } else {
    image(seq(min(OCN$FD$X),max(OCN$FD$X),OCN$cellsize),
               seq(min(OCN$FD$Y),max(OCN$FD$Y),OCN$cellsize),
               t(Zmat),col=colPalette,xlab=" ",ylab=" ",axes=FALSE,asp=1)
  }
  invisible()
}