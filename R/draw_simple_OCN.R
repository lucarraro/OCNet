
draw_simple_OCN <- function(OCN,
                            A_thr_draw=0.002*OCN$dimX*OCN$dimY*OCN$cellsize^2,
                            plot_title="",
                            river_color="#0066FF",
                            EasyDraw=NULL){
  
  if (is.null(EasyDraw)){
  if (OCN$FD$Nnodes>4e4) {
    EasyDraw=TRUE
  } else {EasyDraw=FALSE}
  }
  
  ## plot final state
  AvailableNodes <- setdiff(1:OCN$FD$Nnodes,OCN$FD$Outlet)
  par(pty="s",bty="n",mar=c(1,1,1,1))
  plot(c(min(OCN$FD$X),max(OCN$FD$X)),c(min(OCN$FD$Y),max(OCN$FD$Y)),
       type="n",main=plot_title,asp=1,axes=FALSE,xlab=" ",ylab=" ")
  #points(OCN$FD$X[OCN$FD$Outlet],OCN$FD$Y[OCN$FD$Outlet],pch=22,col="#000000",bg="#000000")
  
  if (EasyDraw==FALSE){
    for (i in AvailableNodes){
      if (OCN$FD$A[i]<A_thr_draw  & abs(OCN$FD$X[i]-OCN$FD$X[OCN$FD$DownNode[i]]) <= OCN$cellsize & abs(OCN$FD$Y[i]-OCN$FD$Y[OCN$FD$DownNode[i]]) <= OCN$cellsize  ) {
        lines(c(OCN$FD$X[i],OCN$FD$X[OCN$FD$DownNode[i]]),c(OCN$FD$Y[i],OCN$FD$Y[OCN$FD$DownNode[i]]),lwd=0.5,col="#E0E0E0")}
    } 
  } 
  for (i in AvailableNodes){
    if (OCN$FD$A[i]>=A_thr_draw  & abs(OCN$FD$X[i]-OCN$FD$X[OCN$FD$DownNode[i]]) <= OCN$cellsize & abs(OCN$FD$Y[i]-OCN$FD$Y[OCN$FD$DownNode[i]]) <= OCN$cellsize) {
      lines(c(OCN$FD$X[i],OCN$FD$X[OCN$FD$DownNode[i]]),c(OCN$FD$Y[i],OCN$FD$Y[OCN$FD$DownNode[i]]),
            lwd=0.5+4.5*(OCN$FD$A[i]/(OCN$FD$Nnodes*OCN$cellsize^2))^0.5,col=river_color)}}

}