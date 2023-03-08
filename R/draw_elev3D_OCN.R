
draw_elev3D_OCN <- function(OCN,
                            coarseGrain=c(1,1),
                            colPalette=terrain.colors(1000,alpha=1),
                            addColorbar=TRUE,
                            drawRiver=TRUE,
                            thrADraw=0.002*OCN$FD$nNodes*OCN$cellsize^2,
                            riverColor="#00CCFF",
                            theta=-20,
                            phi=30,
                            expand=0.05,
                            shade=0.5,
                            min_lwd=0.5,
                            max_lwd=5,
                            args_imagePlot = list()){

  if (!("Z" %in% names(OCN$FD))){
    stop('Missing fields in OCN. You should run landscape_OCN prior to draw_elev3D_OCN.')
  }
  if ((OCN$dimX %% coarseGrain[1] != 0) || (OCN$dimY %% coarseGrain[2] != 0)){
    stop('coarseGrain[1] must be divisor of dimX; coarseGrain[2] must be divisor of dimY')
  }  
  
  if (length(OCN$xllcorner)==0){xllcorner <- min(OCN$FD$X)[1]} else {xllcorner <- OCN$xllcorner}
  if (length(OCN$yllcorner)==0){yllcorner <- min(OCN$FD$Y)[1]} else {yllcorner <- OCN$yllcorner}
  
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
  
  Xvec <- seq(xllcorner,xllcorner+(OCN$dimX-1)*OCN$cellsize,OCN$cellsize)
  Yvec <- seq(yllcorner,yllcorner+(OCN$dimY-1)*OCN$cellsize,OCN$cellsize)
  
  Z_cg <- matrix(data=0,nrow=OCN$dimY/coarseGrain[2],ncol=OCN$dimX/coarseGrain[1])
  X_cg <- rep(0,OCN$dimX/coarseGrain[1])
  Y_cg <- rep(0,OCN$dimY/coarseGrain[2])
  
  for (i in 1:(OCN$dimX/coarseGrain[1])){
    subX <- ((i-1)*coarseGrain[1]+1):(i*coarseGrain[1])
    X_cg[i] <- mean(Xvec[subX])
    for (j in 1:(OCN$dimY/coarseGrain[2])){
      subY <- ((j-1)*coarseGrain[2]+1):(j*coarseGrain[2])
      Z_cg[j,i] <- mean(Zmat[subY,subX])
      Y_cg[j] <- mean(Yvec[subY])
    }
  }
  Zmat <- Z_cg
  
  #old.par <- par(no.readonly = TRUE)
  #on.exit(par(old.par))
  #par(bty="n")
  # draw 3D elevation map
  #Zmat<-matrix(data=OCN$FD$Z,nrow=OCN$dimY,ncol=OCN$dimX)
  zfacet <- Zmat[-1, -1] + Zmat[-1, -OCN$dimX/coarseGrain[1]] + Zmat[-OCN$dimY/coarseGrain[2], -1] + Zmat[-OCN$dimY/coarseGrain[2], -OCN$dimX/coarseGrain[1]]
  vt <- persp(seq(min(X_cg),max(X_cg),OCN$cellsize*coarseGrain[1]),seq(min(Y_cg),max(Y_cg),OCN$cellsize*coarseGrain[2]),t(Zmat),
              theta=theta,phi=phi,col=colPalette[cut(t(zfacet),1000)],expand=expand,shade=shade,
              border=NA,axes=FALSE)
  
  if (addColorbar==TRUE){
  if(is.null(args_imagePlot$col)){args_imagePlot$col <- colPalette}
  if(is.null(args_imagePlot$zlim)){args_imagePlot$zlim <- c(min(OCN$FD$Z),max(OCN$FD$Z))} 
    args_imagePlot$legend.only <- TRUE
    do.call(imagePlot, args_imagePlot)
  #imagePlot(col=colPalette,legend.only=TRUE,zlim=c(min(OCN$FD$Z),max(OCN$FD$Z)))
  }
  # draw river
  if (drawRiver==TRUE){
    river <- trans3d(OCN$FD$X,OCN$FD$Y,OCN$FD$Z,vt)
    AvailableNodes <- setdiff(1:OCN$FD$nNodes,OCN$FD$outlet)
    for (i in AvailableNodes){
      if (OCN$FD$A[i]>=thrADraw & 
          abs(OCN$FD$X[i]-OCN$FD$X[OCN$FD$downNode[i]]) <= 1.001*OCN$cellsize & 
          abs(OCN$FD$Y[i]-OCN$FD$Y[OCN$FD$downNode[i]]) <= 1.001*OCN$cellsize) {
        lines(c(river[[1]][i],river[[1]][OCN$FD$downNode[i]]),c(river[[2]][i],river[[2]][OCN$FD$downNode[i]]),
              lwd=min_lwd+(max_lwd-min_lwd)*(OCN$FD$A[i]/(OCN$FD$nNodes*OCN$cellsize^2))^0.5,col=riverColor)}
    }
  }
  invisible()
}