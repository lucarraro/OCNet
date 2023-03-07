
draw_elev3Drgl_OCN <- function(OCN,
                               coarseGrain=c(1,1),
                               chooseCM=FALSE,
                               addColorbar=FALSE,
                               drawRiver=FALSE,
                               thrADraw=0.002*OCN$FD$nNodes*OCN$cellsize^2,
                               riverColor="#00CCFF",
                               min_lwd=1,
                               max_lwd=8,
                               ...){
  #aspect=c(1,1,0.1),
  #ColPalette=terrain.colors(1000,alpha=1),
  
  
  # give default values to unspecified arguments
  args.def <- list(aspect=c(OCN$dimX/sqrt(OCN$dimX*OCN$dimY),OCN$dimY/sqrt(OCN$dimX*OCN$dimY),0.1),axes=FALSE,xlab="",ylab="",zlab="")
  inargs <- list(...)
  args.def[names(inargs)] <- inargs
  
  if (!("Z" %in% names(OCN$FD))){
    stop('Missing fields in OCN. You should run landscape_OCN prior to draw_elev3D_OCN.')
  }
  
  if (!(chooseCM %in% 1:length(OCN$CM$A)) && !is.logical(chooseCM)) {
    stop('Invalid choice for chooseCM.')
  }
  
  if ((OCN$dimX %% coarseGrain[1] != 0) || (OCN$dimY %% coarseGrain[2] != 0)){
    stop('coarseGrain[1] must be divisor of dimX; coarseGrain[2] must be divisor of dimY')
  }   
  
  if ( (!chooseCM || OCN$CM$A[chooseCM] == OCN$nNodes*OCN$cellsize^2) && OCN$FD$nNodes==OCN$dimX*OCN$dimY){ 
    
    if (is.null(OCN$xllcorner)){xllcorner <- min(OCN$FD$X)[1]} else {xllcorner <- OCN$xllcorner}
    if (is.null(OCN$yllcorner)){yllcorner <- min(OCN$FD$Y)[1]} else {yllcorner <- OCN$yllcorner}
    
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
        Z_cg[j,i] <- mean(Zmat[subY,subX],na.rm=T)
        Y_cg[j] <- mean(Yvec[subY])
      }
    }
    
    par3d(windowRect = c(0, 30, 1000, 1000))
    Zmat <- Zmat + 0.005*mean(Zmat,na.rm=T)
    zlim <- range(Zmat, na.rm=T)
    
    zlim[1] <- floor(zlim[1]- 0.005*mean(Zmat,na.rm=T)); zlim[2] <- ceiling(zlim[2]) 
    zlen <- zlim[2] - zlim[1] + 1
    colorlut <- terrain.colors(zlen) # height color lookup table
    col <- colorlut[ t(Z_cg) - zlim[1] + 1 ] # assign colors to heights for each point
    
    #persp3d(X_cg,Y_cg,t(Z_cg), color = col, zlim=zlim,axes=FALSE,xlab="",ylab="",zlab="",aspect=arguments$aspect,...) # 
    
    do.call(persp3d,c(list(x=X_cg,y=Y_cg,z=t(Z_cg),color = col, zlim=zlim), args.def))
    polygon3d(c(X_cg[1],X_cg[1],X_cg[length(X_cg)],X_cg[length(X_cg)]), c(Y_cg[1],Y_cg[length(Y_cg)],Y_cg[length(Y_cg)],Y_cg[1]), 
              c(rep(zlim[1],4)),coords=c(1,2),col="#997070")
    
    suppressWarnings(try(polygon3d(c(X_cg[1],rep(X_cg[1],length(Y_cg)),X_cg[1]), c(Y_cg[1],Y_cg,Y_cg[length(Y_cg)]), 
                                   c(zlim[1],Z_cg[,1],zlim[1]),coords=2:3,col="#997070"),silent=TRUE))
    suppressWarnings(try(polygon3d(c(X_cg[1],X_cg,X_cg[length(X_cg)]), c(Y_cg[1],rep(Y_cg[1],length(X_cg)),Y_cg[1]), 
                                   c(zlim[1],Z_cg[1,],zlim[1]),coords=c(1,3),col="#997070"),silent=TRUE))
    suppressWarnings(try(polygon3d(c(X_cg[length(X_cg)],rep(X_cg[length(X_cg)],length(Y_cg)),X_cg[length(X_cg)]), c(Y_cg[1],Y_cg,Y_cg[length(Y_cg)]), 
                                   c(zlim[1],Z_cg[,length(X_cg)],zlim[1]),coords=2:3,col="#997070"),silent=TRUE))
    suppressWarnings(try(polygon3d(c(X_cg[1],X_cg,X_cg[length(X_cg)]), c(Y_cg[length(Y_cg)],rep(Y_cg[length(Y_cg)],length(X_cg)),Y_cg[length(Y_cg)]), 
                                   c(zlim[1],Z_cg[length(Y_cg),],zlim[1]),coords=c(1,3),col="#997070"),silent=TRUE))
    
    
    offset <- 0.1*(zlim[2]-zlim[1])
    # draw river
    if (drawRiver==TRUE){
      AvailableNodes <- setdiff(1:OCN$FD$nNodes,OCN$FD$outlet)
      for (i in AvailableNodes){
        if (OCN$FD$A[i]>=thrADraw & 
            abs(OCN$FD$X[i]-OCN$FD$X[OCN$FD$downNode[i]]) <= 1.001*OCN$cellsize & 
            abs(OCN$FD$Y[i]-OCN$FD$Y[OCN$FD$downNode[i]]) <= 1.001*OCN$cellsize) {
          lines3d(c(OCN$FD$X[i],OCN$FD$X[OCN$FD$downNode[i]]),c(OCN$FD$Y[i],OCN$FD$Y[OCN$FD$downNode[i]]),
                  offset + c(OCN$FD$Z[i],OCN$FD$Z[OCN$FD$downNode[i]]),
                  lwd=min_lwd+(max_lwd-min_lwd)*(OCN$FD$A[i]/(OCN$FD$nNodes*OCN$cellsize^2))^0.5,col=riverColor)}
      }
    }
    # add colorbar
    if (addColorbar==TRUE){  
      bgplot3d(suppressWarnings(imagePlot(legend.only=TRUE, zlim=zlim,col=colorlut)))
    }
    
  } else {
    
    if (is.logical(chooseCM)){
      chooseCM <- which(OCN$CM$A==max(OCN$CM$A))
    }
    
    if (OCN$periodicBoundaries==FALSE){
      
      mask <- which(OCN$FD$toCM!=chooseCM)
      
      if (is.null(OCN$xllcorner)){xllcorner <- min(OCN$FD$X)[1]} else {xllcorner <- OCN$xllcorner}
      if (is.null(OCN$yllcorner)){yllcorner <- min(OCN$FD$Y)[1]} else {yllcorner <- OCN$yllcorner}
      
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
      
      Zmat2 <- Zmat
      Zmat2[mask] <- NaN
      
      par3d(windowRect = c(0, 30, 1000, 1000))
      zlim <- range(Zmat,na.rm = T)
      zlim[1] <- floor(zlim[1]); zlim[2] <- ceiling(zlim[2]) 
      zlen <- zlim[2] - zlim[1] + 1
      colorlut <- terrain.colors(zlen) # height color lookup table
      col <- colorlut[ t(Zmat) - zlim[1] + 1 ] # assign colors to heights for each point
      
      persp3d(Xvec,Yvec,t(Zmat2), color = col, zlim=zlim,axes=FALSE,xlab="",ylab="",zlab="",aspect=c(1,1,0.1)) #
      
      if (drawRiver==TRUE){
        offset <- 0.1*(zlim[2]-zlim[1])
        AvailableNodes <- setdiff(1:OCN$FD$nNodes,OCN$FD$outlet)
        for (i in AvailableNodes){
          if (OCN$FD$A[i]>=thrADraw & OCN$FD$toCM[i] == chooseCM &
              abs(OCN$FD$X[i]-OCN$FD$X[OCN$FD$downNode[i]]) <= 1.001*OCN$cellsize & 
              abs(OCN$FD$Y[i]-OCN$FD$Y[OCN$FD$downNode[i]]) <= 1.001*OCN$cellsize) {
            lines3d(c(OCN$FD$X[i],OCN$FD$X[OCN$FD$downNode[i]]),c(OCN$FD$Y[i],OCN$FD$Y[OCN$FD$downNode[i]]),
                    offset + c(OCN$FD$Z[i],OCN$FD$Z[OCN$FD$downNode[i]]),
                    lwd=0.5+7*(OCN$FD$A[i]/(OCN$FD$nNodes*OCN$cellsize^2))^0.5,col="#00CCFF")}
        }
      }
      
      index_border <- NULL
      for (k in 1:length(OCN$CM$XContour[[chooseCM]][[1]])){
        
        x <- OCN$CM$XContour[[chooseCM]][[1]][k] 
        y <- OCN$CM$YContour[[chooseCM]][[1]][k]
        
        ind <- which(abs(OCN$FD$X-x) < OCN$cellsize/2 & abs(OCN$FD$Y-y) < OCN$cellsize/2)
        index_border <- c(index_border,ind)
      }
      index_border <- unique(index_border)
      
      for (i in 1:(length(index_border)-1)) {
        quads3d(c(OCN$FD$X[index_border[i]],OCN$FD$X[index_border[i+1]],OCN$FD$X[index_border[i+1]],OCN$FD$X[index_border[i]]),
                c(OCN$FD$Y[index_border[i]],OCN$FD$Y[index_border[i+1]],OCN$FD$Y[index_border[i+1]],OCN$FD$Y[index_border[i]]),
                c(0,0,OCN$FD$Z[index_border[i+1]],OCN$FD$Z[index_border[i]]),col="#997070")
      }
      quads3d(c(OCN$FD$X[index_border[length(index_border)]],OCN$FD$X[index_border[1]],OCN$FD$X[index_border[1]],OCN$FD$X[index_border[length(index_border)]]),
              c(OCN$FD$Y[index_border[length(index_border)]],OCN$FD$Y[index_border[1]],OCN$FD$Y[index_border[1]],OCN$FD$Y[index_border[length(index_border)]]),
              c(0,0,OCN$FD$Z[index_border[1]],OCN$FD$Z[index_border[length(index_border)]]),col="#997070")
      
      
      lines3d(OCN$FD$X[index_border],OCN$FD$Y[index_border],OCN$FD$Z[index_border],col="black")
      polygon3d(OCN$FD$X[index_border],OCN$FD$Y[index_border],numeric(length(index_border)),col="#997070")
      lines3d(c(OCN$FD$X[index_border],OCN$FD$X[index_border[1]]),c(OCN$FD$Y[index_border],OCN$FD$Y[index_border[1]]),numeric(length(index_border)+1),col="black")
      
      if (addColorbar==TRUE){  
        bgplot3d(suppressWarnings(imagePlot(legend.only=TRUE, zlim=zlim,col=colorlut)))
      }
      
      
    } else { # use "draw" coordinates if periodicBoundaries = TRUE
      
      subset_X <- OCN$FD$XDraw[OCN$FD$toCM==chooseCM]
      subset_Y <- OCN$FD$YDraw[OCN$FD$toCM==chooseCM]
      X_mesh <- seq(min(subset_X)-OCN$cellsize,max(subset_X)+OCN$cellsize,OCN$cellsize)
      Y_mesh <- seq(min(subset_Y)-OCN$cellsize,max(subset_Y)+OCN$cellsize,OCN$cellsize)
      Z_mesh <- matrix(data=NaN,nrow=length(Y_mesh),ncol=length(X_mesh))
      for (i in which(OCN$FD$toCM==chooseCM)){
        ind_x <- which(X_mesh==OCN$FD$XDraw[i])
        ind_y <- which(Y_mesh==OCN$FD$YDraw[i])
        Z_mesh[ind_y,ind_x] <- OCN$FD$Z[i]
      }
      
      par3d(windowRect = c(0, 30, 1000, 1000))
      zlim <- range(OCN$FD$Z[OCN$FD$toCM==chooseCM])
      zlim[1] <- floor(zlim[1]); zlim[2] <- ceiling(zlim[2]) 
      zlen <- zlim[2] - zlim[1] + 1
      colorlut <- terrain.colors(zlen) # height color lookup table
      col <- colorlut[ t(Z_mesh) - zlim[1] + 1 ] # assign colors to heights for each point
      persp3d(X_mesh,Y_mesh,t(Z_mesh), color = col, zlim=zlim,axes=FALSE,xlab="",ylab="",zlab="",aspect=c(1,1,0.1)) #
      
      index_border <- NULL
      for (k in 1:length(OCN$CM$XContourDraw[[chooseCM]][[1]])){
        
        x <- OCN$CM$XContourDraw[[chooseCM]][[1]][k] 
        y <- OCN$CM$YContourDraw[[chooseCM]][[1]][k]
        
        ind <- which(abs(OCN$FD$XDraw-x) < OCN$cellsize/2 & abs(OCN$FD$YDraw-y) < OCN$cellsize/2)
        index_border <- c(index_border,ind)
      }
      index_border <- unique(index_border)
      
      for (i in 1:(length(index_border)-1)) {
        quads3d(c(OCN$FD$XDraw[index_border[i]],OCN$FD$XDraw[index_border[i+1]],OCN$FD$XDraw[index_border[i+1]],OCN$FD$XDraw[index_border[i]]),
                c(OCN$FD$YDraw[index_border[i]],OCN$FD$YDraw[index_border[i+1]],OCN$FD$YDraw[index_border[i+1]],OCN$FD$YDraw[index_border[i]]),
                c(0,0,OCN$FD$Z[index_border[i+1]],OCN$FD$Z[index_border[i]]),col="#997070")
      }
      quads3d(c(OCN$FD$XDraw[index_border[length(index_border)]],OCN$FD$XDraw[index_border[1]],OCN$FD$XDraw[index_border[1]],OCN$FD$XDraw[index_border[length(index_border)]]),
              c(OCN$FD$YDraw[index_border[length(index_border)]],OCN$FD$YDraw[index_border[1]],OCN$FD$YDraw[index_border[1]],OCN$FD$YDraw[index_border[length(index_border)]]),
              c(0,0,OCN$FD$Z[index_border[1]],OCN$FD$Z[index_border[length(index_border)]]),col="#997070")
      
      lines3d(OCN$FD$XDraw[index_border],OCN$FD$YDraw[index_border],OCN$FD$Z[index_border],col="black")
      polygon3d(OCN$FD$XDraw[index_border],OCN$FD$YDraw[index_border],numeric(length(index_border)),col="#997070")
      lines3d(c(OCN$FD$XDraw[index_border],OCN$FD$XDraw[index_border[1]]),c(OCN$FD$YDraw[index_border],OCN$FD$YDraw[index_border[1]]),numeric(length(index_border)+1),col="black")
      
      if (drawRiver==TRUE){
        offset <- 0.1*(zlim[2]-zlim[1])
        AvailableNodes <- setdiff(1:OCN$FD$nNodes,OCN$FD$outlet)
        for (i in AvailableNodes){
          if (OCN$FD$A[i]>=thrADraw & OCN$FD$toCM[i] == chooseCM &
              abs(OCN$FD$XDraw[i]-OCN$FD$XDraw[OCN$FD$downNode[i]]) <= OCN$cellsize & 
              abs(OCN$FD$YDraw[i]-OCN$FD$YDraw[OCN$FD$downNode[i]]) <= OCN$cellsize) {
            lines3d(c(OCN$FD$XDraw[i],OCN$FD$XDraw[OCN$FD$downNode[i]]),c(OCN$FD$YDraw[i],OCN$FD$YDraw[OCN$FD$downNode[i]]),
                    offset + c(OCN$FD$Z[i],OCN$FD$Z[OCN$FD$downNode[i]]),
                    lwd=0.5+7*(OCN$FD$A[i]/(OCN$FD$nNodes*OCN$cellsize^2))^0.5,col="#00CCFF")}
        }
      }
      
    }
  }
}