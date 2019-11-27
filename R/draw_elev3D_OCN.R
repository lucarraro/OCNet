
draw_elev3D_OCN <- function(OCN,
                            coarse_grain=c(1,1),
                            ColPalette=terrain.colors(1000,alpha=1),
                            add_colorbar=TRUE,
                            draw_river=TRUE,
                            A_thr_draw=0.002*OCN$dimX*OCN$dimY*OCN$cellsize^2,
                            river_color="#00CCFF",
                            theta=-20,
                            phi=30,
                            expand=0.05,
                            shade=0.5){

  if (!("Z" %in% names(OCN$FD))){
    stop('Missing fields in OCN. You should run landscape_OCN prior to draw_elev3D_OCN.')
  }
  
  Zmat <- matrix(data=OCN$FD$Z,nrow=OCN$dimY,ncol=OCN$dimX)
  Xvec <- seq(min(OCN$FD$X),max(OCN$FD$X),OCN$cellsize)
  Yvec <- seq(min(OCN$FD$Y),max(OCN$FD$Y),OCN$cellsize)
  
  Z_cg <- matrix(data=0,nrow=OCN$dimY/coarse_grain[2],ncol=OCN$dimX/coarse_grain[1])
  X_cg <- rep(0,OCN$dimX/coarse_grain[1])
  Y_cg <- rep(0,OCN$dimY/coarse_grain[2])
  
  for (i in 1:(OCN$dimX/coarse_grain[1])){
    subX <- ((i-1)*coarse_grain[1]+1):(i*coarse_grain[1])
    X_cg[i] <- mean(Xvec[subX])
    for (j in 1:(OCN$dimY/coarse_grain[2])){
      subY <- ((j-1)*coarse_grain[2]+1):(j*coarse_grain[2])
      Z_cg[j,i] <- mean(Zmat[subY,subX])
      Y_cg[j] <- mean(Yvec[subY])
    }
  }
  Zmat <- Z_cg
  
  par(bty="n",mar=c(1,1,1,1))
  # draw 3D elevation map
  #Zmat<-matrix(data=OCN$FD$Z,nrow=OCN$dimY,ncol=OCN$dimX)
  zfacet <- Zmat[-1, -1] + Zmat[-1, -OCN$dimX/coarse_grain(1)] + Zmat[-OCN$dimY/coarse_grain(2), -1] + Zmat[-OCN$dimY/coarse_grain(2), -OCN$dimX/coarse_grain(1)]
  vt <- persp(seq(min(X_cg),max(X_cg),OCN$cellsize*coarse_grain(1)),seq(min(Y_cg),max(Y_cg),OCN$cellsize*coarse_grain(2)),t(Zmat),
              theta=theta,phi=phi,col=ColPalette[cut(t(zfacet),1000)],expand=expand,shade=shade,
              border=NA,axes=FALSE)
  
  if (add_colorbar==TRUE){
  image.plot(col=ColPalette,legend.only=TRUE,zlim=c(min(OCN$FD$Z),max(OCN$FD$Z)))
  }
  # draw river
  if (draw_river==TRUE){
    river <- trans3d(OCN$FD$X,OCN$FD$Y,OCN$FD$Z,vt)
    AvailableNodes <- setdiff(1:OCN$FD$Nnodes,OCN$FD$Outlet)
    for (i in AvailableNodes){
      if (OCN$FD$A[i]>=A_thr_draw & 
          abs(OCN$FD$X[i]-OCN$FD$X[OCN$FD$DownNode[i]]) <= OCN$cellsize & 
          abs(OCN$FD$Y[i]-OCN$FD$Y[OCN$FD$DownNode[i]]) <= OCN$cellsize) {
        lines(c(river[[1]][i],river[[1]][OCN$FD$DownNode[i]]),c(river[[2]][i],river[[2]][OCN$FD$DownNode[i]]),
              lwd=0.5+4.5*(OCN$FD$A[i]/(OCN$FD$Nnodes*OCN$cellsize^2))^0.5,col=river_color)}
    }
  }
  
}