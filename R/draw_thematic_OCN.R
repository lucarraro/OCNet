
draw_thematic_OCN <- function(theme,OCN,
                              minval=min(theme[!(is.nan(theme))]),
                              maxval=max(theme[!(is.nan(theme))]),
                              plot_title="",
                              ColPalette=colorRampPalette(c("yellow","red","black")),
                              ColLevels=1000,
                              DiscreteLevels=FALSE,
                              ExactDraw=FALSE,
                              DrawNodes=FALSE,
                              Node_cex=2,
                              Node_pch=21,
                              NaNColor="#0099FF",
                              BackgroundColor=NULL,
                              add_colorbar=TRUE){
  
  if (!("RN" %in% names(OCN))){
    stop('Missing fields in OCN. You should run aggregate_OCN prior to draw_thematic_OCN.')
  }
  
  if (minval==maxval){maxval=minval+1}
  
  if (DiscreteLevels==TRUE) {ColLevels <- length(unique(theme))}
  
  if (typeof(ColPalette)=="closure") {
    ColPalette <- ColPalette(ColLevels)
  } else if (typeof(ColPalette)=="character") {
    ColPalette <- ColPalette[1:ColLevels] }
  
  tmp <- "TMP"
  if (length(theme)==OCN$RN$Nnodes && (length(theme)==OCN$AG$Nnodes)){
    while ((tmp != "RN") && (tmp != "AG"))
      tmp <- readline(prompt="theme can be interpreted as a vector both at the RN and AG levels. Choose desired level by typing RN or AG: ")
    if (tmp == "RN"){
      byRN = TRUE
    } else if (tmp == "AG"){
      byRN = FALSE
    } else {
      print('Wrong input!')
    }
  } 
  
  if (length(theme)==OCN$RN$Nnodes){
    byRN = TRUE
  } else if (length(theme)==OCN$AG$Nnodes){
    byRN = FALSE
  } else {
    stop('theme has invalid length')
  }
  
  if (length(Node_cex)>1 && length(Node_cex) != length(theme)){
    stop('Node_cex has invalid length')
  }
  
  if (length(Node_pch)>1 && length(Node_pch) != length(theme)){
    stop('Node_pch has invalid length')
  }
  
  if (ExactDraw==TRUE){
    X <- OCN$FD$X_draw
    Y <- OCN$FD$Y_draw
    Xc <- OCN$CM$X_contour_draw
    Yc <- OCN$CM$Y_contour_draw
  } else {
    X <- OCN$FD$X
    Y <- OCN$FD$Y
    Xc <- OCN$CM$X_contour
    Yc <- OCN$CM$Y_contour
  }
  
  if (length(Node_cex)==1){
    cex_vec <- Node_cex*rep(1,length(theme))
  } else {cex_vec <- Node_cex}
  
  if (length(Node_pch)==1){
    pch_vec <- Node_pch*rep(1,length(theme))
  } else {pch_vec <- Node_pch}
  
  
  AvailableNodes <- setdiff(1:OCN$FD$Nnodes,OCN$FD$Outlet)
  par(bty="n",mar=c(1,1,1,1))
  plot(c(min(X),max(X)),c(min(Y),max(Y)),
       type="n",xlab=" ",ylab=" ",main=plot_title,axes=FALSE,asp=1)
  
  if (!is.null(BackgroundColor)){
    if ((OCN$N_outlet > 1) && (length(BackgroundColor)==1) ){
      BackgroundColor=rep(BackgroundColor,OCN$N_outlet)}
    for (i in 1:OCN$N_outlet){
      for (j in 1:length(Xc[[i]])){
    polygon(Xc[[i]][[j]],Yc[[i]][[j]],col=BackgroundColor[i],lty=0)
      }}
  }
  #for (j in 1:OCN$N_outlet){
  #  for (k in 1:length(Xc[[j]])){
  #    lines(Xc[[j]][[k]],Yc[[j]][[k]],lwd=2)
  #  }
  #}
  #points(X[OCN$FD$Outlet],Y[OCN$FD$Outlet],pch=22,col="#000000",bg="#000000")
  
  for (i in AvailableNodes){
    rn <- OCN$FD$to_RN[i]
    reach <- OCN$RN$to_AG[OCN$FD$to_RN[i]]
    if (OCN$FD$A[i]>=OCN$A_thr & 
        abs(X[i]-X[OCN$FD$DownNode[i]]) <= OCN$cellsize & 
        abs(Y[i]-Y[OCN$FD$DownNode[i]]) <= OCN$cellsize  ) {
      if (is.nan(theme[reach])==TRUE | is.na(theme[reach])==TRUE) {
        hexcolor <- NaNColor
      } else {
        if (byRN==TRUE){
          colvalue <- 1+round((ColLevels-1)*max(0,min(1,(theme[rn]-minval)/(maxval-minval))))
        } else {colvalue <- 1+round((ColLevels-1)*max(0,min(1,(theme[reach]-minval)/(maxval-minval))))}
        
        hexcolor <- ColPalette[colvalue]
      }
      if (DrawNodes==TRUE){
        lines(c(X[i],X[OCN$FD$DownNode[i]]),c(Y[i],Y[OCN$FD$DownNode[i]]),
              lwd=0.5+4.5*(OCN$FD$A[i]/(OCN$FD$Nnodes*OCN$cellsize^2))^0.5,col=NaNColor)
      } else {
        lines(c(X[i],X[OCN$FD$DownNode[i]]),c(Y[i],Y[OCN$FD$DownNode[i]]),
              lwd=0.5+4.5*(OCN$FD$A[i]/(OCN$FD$Nnodes*OCN$cellsize^2))^0.5,col=hexcolor)}
    }
  }
  
  if (DrawNodes==TRUE && ExactDraw==FALSE){
    if (byRN==TRUE){
      Nnodes <- OCN$RN$Nnodes
    } else {Nnodes <- OCN$AG$Nnodes}
    for (i in 1:Nnodes){
      colvalue <- 1+round((ColLevels-1)*max(0,min(1,(theme[i]-minval)/(maxval-minval))))
      hexcolor <- ColPalette[colvalue]
      if (byRN==TRUE){
        points(OCN$RN$X[i],OCN$RN$Y[i],bg=hexcolor,pch=pch_vec[i],cex=cex_vec[i])
      } else {
        points(OCN$AG$X[i],OCN$AG$Y[i],bg=hexcolor,pch=pch_vec[i],cex=cex_vec[i])  
      }
    }
  }
  
  if (DrawNodes==TRUE && ExactDraw==TRUE){
    if (byRN==TRUE){
      Nnodes <- OCN$RN$Nnodes
    } else {Nnodes <- OCN$AG$Nnodes}
    for (i in 1:Nnodes){
      colvalue <- 1+round((ColLevels-1)*max(0,min(1,(theme[i]-minval)/(maxval-minval))))
      hexcolor <- ColPalette[colvalue]
      if (byRN==TRUE){
        node <- which(OCN$FD$X==OCN$RN$X[i] & OCN$FD$Y==OCN$RN$Y[i])
        points(X[node],Y[node],bg=hexcolor,pch=pch_vec[i],cex=cex_vec[i])
      } else {
        node <- which(OCN$FD$X==OCN$AG$X[i] & OCN$FD$Y==OCN$AG$Y[i])
        points(X[node],Y[node],bg=hexcolor,pch=pch_vec[i],cex=cex_vec[i])  
      }
    }
  }
  
  
  if (add_colorbar) {image.plot(col=ColPalette,legend.only=TRUE,zlim=c(minval,maxval))}
  
}