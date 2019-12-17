
draw_elev2D_OCN <- function(OCN,
                              colPalette=terrain.colors(1000,alpha=1)){
  
  if (!("Z" %in% names(OCN$FD))){
    stop('Missing fields in OCN. You should run landscape_OCN prior to draw_elev2D_OCN.')
  }
  
  # plot elevation map
  #old.par <- par(no.readonly =TRUE) 
  #on.exit(par(old.par))
  #par(bty="n")
  Zmat<-matrix(data=OCN$FD$Z,nrow=OCN$dimY,ncol=OCN$dimX)
  image.plot(seq(min(OCN$FD$X),max(OCN$FD$X),OCN$cellsize),
             seq(min(OCN$FD$Y),max(OCN$FD$Y),OCN$cellsize),
             t(Zmat),col=colPalette,xlab=" ",ylab=" ",axes=FALSE,asp=1)
  invisible()
}