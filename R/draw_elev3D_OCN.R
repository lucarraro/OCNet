
draw_elev3D_OCN <- function(OCN,
                            ColPalette=terrain.colors(1000,alpha=1),
                            add_colorbar=TRUE,
                            draw_river=TRUE,
                            A_thr_draw=0.002*OCN$dimX*OCN$dimY*OCN$cellsize^2,
                            river_color="#00CCFF"){

  if (!("Z" %in% names(OCN$FD))){
    stop('Missing fields in OCN. You should run landscape_OCN prior to draw_elev3D_OCN.')
  }
  
  par(bty="n",mar=c(1,1,1,1))
  # draw 3D elevation map
  Zmat<-matrix(data=OCN$FD$Z,nrow=OCN$dimY,ncol=OCN$dimX)
  zfacet <- Zmat[-1, -1] + Zmat[-1, -OCN$dimX] + Zmat[-OCN$dimY, -1] + Zmat[-OCN$dimY, -OCN$dimX]
  vt <- persp(seq(min(OCN$FD$X),max(OCN$FD$X),OCN$cellsize),seq(min(OCN$FD$Y),max(OCN$FD$Y),OCN$cellsize),t(Zmat),
              theta=-20,phi=30,col=ColPalette[cut(t(zfacet),1000)],expand=0.05,shade=0.5,
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