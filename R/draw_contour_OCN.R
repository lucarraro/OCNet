
draw_contour_OCN <- function(OCN,
                             thrA.draw=0.002*OCN$dimX*OCN$dimY*OCN$cellsize^2,
                             exactDraw=TRUE,
                             drawContours=TRUE,
                             colPalRiver=NULL,
                             colPalCont="#000000",
                             drawOutlets=0,
                             pch=15,
                             colPalOut="#000000"){
  
  if (!("X_draw" %in% names(OCN$FD))){
    stop('Missing fields in OCN. You should run landscape_OCN prior to draw_contour_OCN.')
  }
  
  if (OCN$FD$Nnodes>4e4) {
    EasyDraw=TRUE
  } else {EasyDraw=FALSE}
  
  
  if (exactDraw==TRUE){
    X_draw <- OCN$FD$X_draw
    Y_draw <- OCN$FD$Y_draw
    Xc <- OCN$CM$X_contour_draw
    Yc <- OCN$CM$Y_contour_draw
  } else {
    X_draw <- OCN$FD$X
    Y_draw <- OCN$FD$Y
    Xc <- OCN$CM$X_contour
    Yc <- OCN$CM$Y_contour
  }
  
  # plot network (draw mode)
  
  if (is.null(colPalRiver)){
    if (OCN$N_outlet > 1){
      rnbw <- hcl.colors(OCN$N_outlet,palette="Dark 3")
      rnbw <- c(rnbw[(round(OCN$N_outlet/2)+1):OCN$N_outlet],rnbw[1:(round(OCN$N_outlet/2)+1)])
    } else {rnbw <- hcl.colors(3,palette="Dark 3")
    rnbw <- rnbw[3]}
    colPalRiver <- rnbw
    
  } else if (typeof(colPalRiver)=="closure") {
    colPalRiver <- colPalRiver(OCN$N_outlet)
  } else if (typeof(colPalRiver)=="character") {
    if (length(colPalRiver)==1){
      colPalRiver <- rep(colPalRiver,OCN$N_outlet)
    }
  }
  
  if (colPalCont==0){
    colPalCont=colPalRiver
  } else if (typeof(colPalCont)=="closure") {
    colPalCont <- colPalCont(OCN$N_outlet)
  } else if (typeof(colPalCont)=="character") {
    if (length(colPalCont)==1){
      colPalCont <- rep(colPalCont,OCN$N_outlet)
    }
  }
  
  if (colPalOut==0) {
    colPalOut=colPalRiver
  } else if (typeof(colPalOut)=="closure") {
    colPalOut <- colPalOut(OCN$N_outlet)
  } else if (typeof(colPalOut)=="character") {
    if (length(colPalOut)==1){
      colPalOut <- rep(colPalOut,OCN$N_outlet)
    }
  }
  
  AvailableNodes <- setdiff(1:OCN$FD$Nnodes,OCN$FD$Outlet)
  par(bty="n",mar=c(1,1,1,1))
  plot(c(min(X_draw),max(X_draw)),c(min(Y_draw),max(Y_draw)),type="n",xlab=" ",ylab=" ",axes=FALSE,asp=1)
  
  if (drawOutlets==1) {
    for (i in 1:OCN$N_outlet){
    points(X_draw[OCN$FD$Outlet[i]],Y_draw[OCN$FD$Outlet[i]],pch=pch,col=colPalOut[i])
    }
  }
  
  if (EasyDraw==FALSE){
    for (i in AvailableNodes){
      if (OCN$FD$A[i]<thrA.draw & abs(X_draw[i]-X_draw[OCN$FD$DownNode[i]]) <= OCN$cellsize & abs(Y_draw[i]-Y_draw[OCN$FD$DownNode[i]]) <= OCN$cellsize) {
        lines(c(X_draw[i],X_draw[OCN$FD$DownNode[i]]),c(Y_draw[i],Y_draw[OCN$FD$DownNode[i]]),lwd=0.5,col="#E0E0E0")} 
    }
  }
  for (i in AvailableNodes){
    if (OCN$FD$A[i]>=thrA.draw & 
        abs(X_draw[i]-X_draw[OCN$FD$DownNode[i]]) <= OCN$cellsize & 
        abs(Y_draw[i]-Y_draw[OCN$FD$DownNode[i]]) <= OCN$cellsize  ) {
      lines(c(X_draw[i],X_draw[OCN$FD$DownNode[i]]),c(Y_draw[i],Y_draw[OCN$FD$DownNode[i]]),
            lwd=0.5+4.5*(OCN$FD$A[i]/(OCN$FD$Nnodes*OCN$cellsize^2))^0.5,col=colPalRiver[OCN$FD$to_CM[i]])}
  }
  if (drawContours){
    for (j in 1:OCN$N_outlet){
      for (k in 1:length(Xc[[j]])){
        lines(Xc[[j]][[k]],Yc[[j]][[k]],lwd=2,col=colPalCont[j])
      }
    }
  }
  if (drawOutlets==2) {
    for (i in 1:OCN$N_outlet){
      points(X_draw[OCN$FD$Outlet[i]],Y_draw[OCN$FD$Outlet[i]],pch=pch,col=colPalOut[i])
    }
  }
  
}