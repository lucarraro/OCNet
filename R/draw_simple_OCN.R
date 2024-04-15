
draw_simple_OCN <- function(OCN,
                            thrADraw=0.002*OCN$FD$nNodes*OCN$cellsize^2,
                            riverColor="#0066FF",
                            easyDraw=NULL,
                            min_lwd=0.5,
                            max_lwd=5,
                            add=FALSE){
  
  if (is.null(easyDraw)){
  if (OCN$FD$nNodes>4e4) {
    easyDraw=TRUE
  } else {easyDraw=FALSE}
  }
  
  if (is.null(dev.list()) & add==TRUE){
    add <- FALSE
    warning("'add' will be ignored as there is no existing plot")
  }
  
  ## plot final state
  AvailableNodes <- setdiff(1:OCN$FD$nNodes,OCN$FD$outlet)
  #old.par <- par(no.readonly = TRUE)
  #on.exit(par(old.par))
  #par(bty="n")
  if (!add){
  plot(c(min(OCN$FD$X),max(OCN$FD$X)),c(min(OCN$FD$Y),max(OCN$FD$Y)),
       type="n",asp=1,axes=FALSE,xlab=" ",ylab=" ")}
  #points(OCN$FD$X[OCN$FD$outlet],OCN$FD$Y[OCN$FD$outlet],pch=22,col="#000000",bg="#000000")
  
  if (easyDraw==FALSE){
    for (i in AvailableNodes){
      if (OCN$FD$A[i]<thrADraw  & abs(OCN$FD$X[i]-OCN$FD$X[OCN$FD$downNode[i]]) <= 1.001*OCN$cellsize & abs(OCN$FD$Y[i]-OCN$FD$Y[OCN$FD$downNode[i]]) <= 1.001*OCN$cellsize  ) {
        lines(c(OCN$FD$X[i],OCN$FD$X[OCN$FD$downNode[i]]),c(OCN$FD$Y[i],OCN$FD$Y[OCN$FD$downNode[i]]),lwd=min_lwd,col="#E0E0E0")}
    } 
  } 
  for (i in AvailableNodes){
    if (OCN$FD$A[i]>=thrADraw  & abs(OCN$FD$X[i]-OCN$FD$X[OCN$FD$downNode[i]]) <= 1.001*OCN$cellsize & abs(OCN$FD$Y[i]-OCN$FD$Y[OCN$FD$downNode[i]]) <= 1.001*OCN$cellsize) {
      lines(c(OCN$FD$X[i],OCN$FD$X[OCN$FD$downNode[i]]),c(OCN$FD$Y[i],OCN$FD$Y[OCN$FD$downNode[i]]),
            lwd=min_lwd+(max_lwd-min_lwd)*(OCN$FD$A[i]/(OCN$FD$nNodes*OCN$cellsize^2))^0.5,col=riverColor)}}
  invisible()
}