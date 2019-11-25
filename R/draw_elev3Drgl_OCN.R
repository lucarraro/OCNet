
draw_elev3Drgl_OCN <- function(OCN,
                               coarse_grain=c(1,1),
                               chooseCM=FALSE,
                               add_colorbar=FALSE,
                               draw_river=FALSE,
                               A_thr_draw=0.002*OCN$dimX*OCN$dimY*OCN$cellsize^2,
                               river_color="#00CCFF",
                               ...){
  #aspect=c(1,1,0.1),
  #ColPalette=terrain.colors(1000,alpha=1),
  
  options(rgl.useNULL = TRUE)
  
  # give default values to unspecified arguments
  args.def <- list(aspect=c(1,1,0.1),axes=FALSE,xlab="",ylab="",zlab="")
  inargs <- list(...)
  args.def[names(inargs)] <- inargs
  
  if (!("Z" %in% names(OCN$FD))){
    stop('Missing fields in OCN. You should run landscape_OCN prior to draw_elev3D_OCN.')
  }
  
  if (!(chooseCM %in% 1:length(OCN$CM$A)) && !is.logical(chooseCM)) {
    stop('Invalid choice for chooseCM.')
  }
  
  if ( !chooseCM || OCN$CM$A[chooseCM] == OCN$dimX*OCN$dimY*OCN$cellsize^2){ 
    
    Zmat <- matrix(data=OCN$FD$Z,nrow=OCN$dimY,ncol=OCN$dimX)
    Xvec <- seq(min(OCN$FD$X),max(OCN$FD$X),OCN$cellsize)
    Yvec <- seq(min(OCN$FD$Y),max(OCN$FD$Y),OCN$cellsize)
    
    Z_cg <- matrix(data=0,nrow=OCN$dimY/coarse_grain[2],ncol=OCN$dimY/coarse_grain[1])
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
    
    par3d(windowRect = c(0, 30, 1000, 1000))
    Zmat <- Zmat + 0.005*mean(Zmat)
    zlim <- range(Zmat)
    
    zlim[1] <- floor(zlim[1]- 0.005*mean(Zmat)); zlim[2] <- ceiling(zlim[2]) 
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
    if (draw_river==TRUE){
      AvailableNodes <- setdiff(1:OCN$FD$Nnodes,OCN$FD$Outlet)
      for (i in AvailableNodes){
        if (OCN$FD$A[i]>=A_thr_draw & 
            abs(OCN$FD$X[i]-OCN$FD$X[OCN$FD$DownNode[i]]) <= OCN$cellsize & 
            abs(OCN$FD$Y[i]-OCN$FD$Y[OCN$FD$DownNode[i]]) <= OCN$cellsize) {
          lines3d(c(OCN$FD$X[i],OCN$FD$X[OCN$FD$DownNode[i]]),c(OCN$FD$Y[i],OCN$FD$Y[OCN$FD$DownNode[i]]),
                  offset + c(OCN$FD$Z[i],OCN$FD$Z[OCN$FD$DownNode[i]]),
                  lwd=1+7*(OCN$FD$A[i]/(OCN$FD$Nnodes*OCN$cellsize^2))^0.5,col=river_color)}
      }
    }
    # add colorbar
    if (add_colorbar==TRUE){  
      bgplot3d(suppressWarnings(image.plot(legend.only=TRUE, zlim=zlim,col=colorlut)))
    }
    
  } else {
    
    if (chooseCM==TRUE && is.logical(chooseCM)){
      chooseCM <- which(OCN$CM$A==max(OCN$CM$A))
    }
    
    if (OCN$PeriodicBoundaries==FALSE){
      
      mask <- which(OCN$FD$to_CM!=chooseCM)
      
      Xvec <- seq(min(OCN$FD$X),max(OCN$FD$X),OCN$cellsize)
      Yvec <- seq(min(OCN$FD$Y),max(OCN$FD$Y),OCN$cellsize)
      
      Zmat <- matrix(data=OCN$FD$Z,nrow=OCN$dimY,ncol=OCN$dimX)
      Zmat2 <- Zmat
      Zmat2[mask] <- NaN
      
      par3d(windowRect = c(0, 30, 1000, 1000))
      zlim <- range(Zmat)
      zlim[1] <- floor(zlim[1]); zlim[2] <- ceiling(zlim[2]) 
      zlen <- zlim[2] - zlim[1] + 1
      colorlut <- terrain.colors(zlen) # height color lookup table
      col <- colorlut[ t(Zmat) - zlim[1] + 1 ] # assign colors to heights for each point
      
      persp3d(Xvec,Yvec,t(Zmat2), color = col, zlim=zlim,axes=FALSE,xlab="",ylab="",zlab="",aspect=c(1,1,0.1)) #
      
      if (draw_river==TRUE){
        offset <- 0.1*(zlim[2]-zlim[1])
        AvailableNodes <- setdiff(1:OCN$FD$Nnodes,OCN$FD$Outlet)
        for (i in AvailableNodes){
          if (OCN$FD$A[i]>=A_thr_draw & OCN$FD$to_CM[i] == chooseCM &
              abs(OCN$FD$X[i]-OCN$FD$X[OCN$FD$DownNode[i]]) <= OCN$cellsize & 
              abs(OCN$FD$Y[i]-OCN$FD$Y[OCN$FD$DownNode[i]]) <= OCN$cellsize) {
            lines3d(c(OCN$FD$X[i],OCN$FD$X[OCN$FD$DownNode[i]]),c(OCN$FD$Y[i],OCN$FD$Y[OCN$FD$DownNode[i]]),
                    offset + c(OCN$FD$Z[i],OCN$FD$Z[OCN$FD$DownNode[i]]),
                    lwd=0.5+7*(OCN$FD$A[i]/(OCN$FD$Nnodes*OCN$cellsize^2))^0.5,col="#00CCFF")}
        }
      }
      
      index_border <- NULL
      for (k in 1:length(OCN$CM$X_contour[[chooseCM]][[1]])){
        
        x <- OCN$CM$X_contour[[chooseCM]][[1]][k] 
        y <- OCN$CM$Y_contour[[chooseCM]][[1]][k]
        
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
      
      if (add_colorbar==TRUE){  
        bgplot3d(suppressWarnings(image.plot(legend.only=TRUE, zlim=zlim,col=colorlut)))
      }
      
      
    } else { # use "draw" coordinates if PeriodicBoundaries = TRUE
      
      subset_X <- OCN$FD$X_draw[OCN$FD$to_CM==chooseCM]
      subset_Y <- OCN$FD$Y_draw[OCN$FD$to_CM==chooseCM]
      X_mesh <- seq(min(subset_X)-OCN$cellsize,max(subset_X)+OCN$cellsize,OCN$cellsize)
      Y_mesh <- seq(min(subset_Y)-OCN$cellsize,max(subset_Y)+OCN$cellsize,OCN$cellsize)
      Z_mesh <- matrix(data=NaN,nrow=length(Y_mesh),ncol=length(X_mesh))
      for (i in which(OCN$FD$to_CM==chooseCM)){
        ind_x <- which(X_mesh==OCN$FD$X_draw[i])
        ind_y <- which(Y_mesh==OCN$FD$Y_draw[i])
        Z_mesh[ind_y,ind_x] <- OCN$FD$Z[i]
      }
      
      par3d(windowRect = c(0, 30, 1000, 1000))
      zlim <- range(OCN$FD$Z[OCN$FD$to_CM==chooseCM])
      zlim[1] <- floor(zlim[1]); zlim[2] <- ceiling(zlim[2]) 
      zlen <- zlim[2] - zlim[1] + 1
      colorlut <- terrain.colors(zlen) # height color lookup table
      col <- colorlut[ t(Z_mesh) - zlim[1] + 1 ] # assign colors to heights for each point
      persp3d(X_mesh,Y_mesh,t(Z_mesh), color = col, zlim=zlim,axes=FALSE,xlab="",ylab="",zlab="",aspect=c(1,1,0.1)) #
      
      index_border <- NULL
      for (k in 1:length(OCN$CM$X_contour_draw[[chooseCM]][[1]])){
        
        x <- OCN$CM$X_contour_draw[[chooseCM]][[1]][k] 
        y <- OCN$CM$Y_contour_draw[[chooseCM]][[1]][k]
        
        ind <- which(abs(OCN$FD$X_draw-x) < OCN$cellsize/2 & abs(OCN$FD$Y_draw-y) < OCN$cellsize/2)
        index_border <- c(index_border,ind)
      }
      index_border <- unique(index_border)
      
      for (i in 1:(length(index_border)-1)) {
        quads3d(c(OCN$FD$X_draw[index_border[i]],OCN$FD$X_draw[index_border[i+1]],OCN$FD$X_draw[index_border[i+1]],OCN$FD$X_draw[index_border[i]]),
                c(OCN$FD$Y_draw[index_border[i]],OCN$FD$Y_draw[index_border[i+1]],OCN$FD$Y_draw[index_border[i+1]],OCN$FD$Y_draw[index_border[i]]),
                c(0,0,OCN$FD$Z[index_border[i+1]],OCN$FD$Z[index_border[i]]),col="#997070")
      }
      quads3d(c(OCN$FD$X_draw[index_border[length(index_border)]],OCN$FD$X_draw[index_border[1]],OCN$FD$X_draw[index_border[1]],OCN$FD$X_draw[index_border[length(index_border)]]),
              c(OCN$FD$Y_draw[index_border[length(index_border)]],OCN$FD$Y_draw[index_border[1]],OCN$FD$Y_draw[index_border[1]],OCN$FD$Y_draw[index_border[length(index_border)]]),
              c(0,0,OCN$FD$Z[index_border[1]],OCN$FD$Z[index_border[length(index_border)]]),col="#997070")
      
      lines3d(OCN$FD$X_draw[index_border],OCN$FD$Y_draw[index_border],OCN$FD$Z[index_border],col="black")
      polygon3d(OCN$FD$X_draw[index_border],OCN$FD$Y_draw[index_border],numeric(length(index_border)),col="#997070")
      lines3d(c(OCN$FD$X_draw[index_border],OCN$FD$X_draw[index_border[1]]),c(OCN$FD$Y_draw[index_border],OCN$FD$Y_draw[index_border[1]]),numeric(length(index_border)+1),col="black")
      
      if (draw_river==TRUE){
        offset <- 0.1*(zlim[2]-zlim[1])
        AvailableNodes <- setdiff(1:OCN$FD$Nnodes,OCN$FD$Outlet)
        for (i in AvailableNodes){
          if (OCN$FD$A[i]>=A_thr_draw & OCN$FD$to_CM[i] == chooseCM &
              abs(OCN$FD$X_draw[i]-OCN$FD$X_draw[OCN$FD$DownNode[i]]) <= OCN$cellsize & 
              abs(OCN$FD$Y_draw[i]-OCN$FD$Y_draw[OCN$FD$DownNode[i]]) <= OCN$cellsize) {
            lines3d(c(OCN$FD$X_draw[i],OCN$FD$X_draw[OCN$FD$DownNode[i]]),c(OCN$FD$Y_draw[i],OCN$FD$Y_draw[OCN$FD$DownNode[i]]),
                    offset + c(OCN$FD$Z[i],OCN$FD$Z[OCN$FD$DownNode[i]]),
                    lwd=0.5+7*(OCN$FD$A[i]/(OCN$FD$Nnodes*OCN$cellsize^2))^0.5,col="#00CCFF")}
        }
      }
      
    }
  }
}