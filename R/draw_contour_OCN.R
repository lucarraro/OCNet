
draw_contour_OCN <- function(OCN,
                             thrADraw=0.002*OCN$dimX*OCN$dimY*OCN$cellsize^2,
                             exactDraw=TRUE,
                             drawContours=TRUE,
                             colPalRiver=NULL,
                             colPalCont="#000000",
                             drawOutlets=0,
                             pch=15,
                             colPalOut="#000000"){
  
  if (!("XDraw" %in% names(OCN$FD))){
    stop('Missing fields in OCN. You should run landscape_OCN prior to draw_contour_OCN.')
  }
  
  if (OCN$FD$nNodes>4e4) {
    EasyDraw=TRUE
  } else {EasyDraw=FALSE}
  
  
  if (exactDraw==TRUE){
    XDraw <- OCN$FD$XDraw
    YDraw <- OCN$FD$YDraw
    Xc <- OCN$CM$XContourDraw
    Yc <- OCN$CM$YContourDraw
  } else {
    XDraw <- OCN$FD$X
    YDraw <- OCN$FD$Y
    Xc <- OCN$CM$XContour
    Yc <- OCN$CM$YContour
  }
  
  # plot network (draw mode)
  
  if (is.null(colPalRiver)){
    if (OCN$nOutlet > 1){
      rnbw <- hcl.colors(OCN$nOutlet,palette="Dark 3")
      rnbw <- c(rnbw[(round(OCN$nOutlet/2)+1):OCN$nOutlet],rnbw[1:(round(OCN$nOutlet/2)+1)])
    } else {rnbw <- hcl.colors(3,palette="Dark 3")
    rnbw <- rnbw[3]}
    colPalRiver <- rnbw
    
  } else if (typeof(colPalRiver)=="closure") {
    colPalRiver <- colPalRiver(OCN$nOutlet)
  } else if (typeof(colPalRiver)=="character") {
    if (length(colPalRiver)==1){
      colPalRiver <- rep(colPalRiver,OCN$nOutlet)
    }
  }
  
  if (colPalCont==0){
    colPalCont=colPalRiver
  } else if (typeof(colPalCont)=="closure") {
    colPalCont <- colPalCont(OCN$nOutlet)
  } else if (typeof(colPalCont)=="character") {
    if (length(colPalCont)==1){
      colPalCont <- rep(colPalCont,OCN$nOutlet)
    }
  }
  
  if (colPalOut==0) {
    colPalOut=colPalRiver
  } else if (typeof(colPalOut)=="closure") {
    colPalOut <- colPalOut(OCN$nOutlet)
  } else if (typeof(colPalOut)=="character") {
    if (length(colPalOut)==1){
      colPalOut <- rep(colPalOut,OCN$nOutlet)
    }
  }
  
  AvailableNodes <- setdiff(1:OCN$FD$nNodes,OCN$FD$outlet)
  old.par <- par(no.readonly =TRUE) 
  on.exit(par(old.par))
  par(bty="n")
  plot(c(min(XDraw),max(XDraw)),c(min(YDraw),max(YDraw)),type="n",xlab=" ",ylab=" ",axes=FALSE,asp=1)
  
  if (drawOutlets==1) {
    for (i in 1:OCN$nOutlet){
    points(XDraw[OCN$FD$outlet[i]],YDraw[OCN$FD$outlet[i]],pch=pch,col=colPalOut[i])
    }
  }
  
  if (EasyDraw==FALSE){
    for (i in AvailableNodes){
      if (OCN$FD$A[i]<thrADraw & abs(XDraw[i]-XDraw[OCN$FD$downNode[i]]) <= OCN$cellsize & abs(YDraw[i]-YDraw[OCN$FD$downNode[i]]) <= OCN$cellsize) {
        lines(c(XDraw[i],XDraw[OCN$FD$downNode[i]]),c(YDraw[i],YDraw[OCN$FD$downNode[i]]),lwd=0.5,col="#E0E0E0")} 
    }
  }
  for (i in AvailableNodes){
    if (OCN$FD$A[i]>=thrADraw & 
        abs(XDraw[i]-XDraw[OCN$FD$downNode[i]]) <= OCN$cellsize & 
        abs(YDraw[i]-YDraw[OCN$FD$downNode[i]]) <= OCN$cellsize  ) {
      lines(c(XDraw[i],XDraw[OCN$FD$downNode[i]]),c(YDraw[i],YDraw[OCN$FD$downNode[i]]),
            lwd=0.5+4.5*(OCN$FD$A[i]/(OCN$FD$nNodes*OCN$cellsize^2))^0.5,col=colPalRiver[OCN$FD$toCM[i]])}
  }
  if (drawContours){
    for (j in 1:OCN$nOutlet){
      for (k in 1:length(Xc[[j]])){
        lines(Xc[[j]][[k]],Yc[[j]][[k]],lwd=2,col=colPalCont[j])
      }
    }
  }
  if (drawOutlets==2) {
    for (i in 1:OCN$nOutlet){
      points(XDraw[OCN$FD$outlet[i]],YDraw[OCN$FD$outlet[i]],pch=pch,col=colPalOut[i])
    }
  }
  invisible()
}