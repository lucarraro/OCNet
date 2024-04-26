
draw_thematic_OCN <- function(OCN,theme=NA*numeric(OCN$AG$nNodes),
                              chooseAggregation=NULL,
                              discreteLevels=FALSE,
                              colLevels=NULL,
                              cutoff=FALSE,
                              colPalette=colorRampPalette(c("yellow","red","black")),
                              exactDraw=FALSE,
                              chooseCM=FALSE,
                              drawNodes=FALSE,
                              nodeType="upstream",
                              nanColor="#00BFFF",
                              riverColor="#00BFFF",
                              backgroundColor="#999999",
                              addLegend=TRUE,
                              min_lwd=0.5,
                              max_lwd=5,
                              add=FALSE,
                              args_imagePlot=list(),
                              args_legend=list(),
                              ...){
  
  dots <- list(...)
  if (is.null(dots$axes)){
    if ("xlim" %in% names(dots) | "ylim" %in% names(dots)) {dots$axes <- TRUE} else {dots$axes <- FALSE}}
  if (is.null(dots$xlab)){dots$xlab <- ""}
  if (is.null(dots$ylab)){dots$ylab <- ""}
  if (is.null(dots$type)){dots$type <- "n"}
  if (is.null(dots$asp)){dots$asp <- 1}
  if (is.null(dots$cex)){dots$cex <- 2}; cex <- dots$cex
  if (is.null(dots$pch)){dots$pch <- 21}; pch <- dots$pch
  
  if (is.null(dev.list()) & add==TRUE){
    add <- FALSE
    warning("'add' will be ignored as there is no existing plot")
  }
  
  # initialization
  if (!("RN" %in% names(OCN))){ # try to swap arguments if in wrong order
    tmp <- OCN
    OCN <- theme
    theme <- tmp
  }
  
  if (all(is.na(theme))) {addLegend <- FALSE }
  
  if (length(OCN$RN$nNodes)==0){stop('Missing fields in OCN. You should run aggregate_OCN prior to draw_thematic_OCN.')}
  
  if (discreteLevels == FALSE) {
    if (length(colLevels)<3){N_colLevels <- 1000} else {N_colLevels <- colLevels[3]}
    
    if (is.null(colLevels)){
      minval <- min(theme[!(is.nan(theme))])
      maxval <- max(theme[!(is.nan(theme))])
      if (is.na(minval) & is.na(maxval)){
        minval <- 0; maxval <- 0;
      }
      colLevels <- c(minval,maxval,N_colLevels)
    }
    minval <- colLevels[1]
    maxval <- colLevels[2]
    if (minval==maxval) {maxval <- minval + 1}
    Breakpoints <- seq(minval,maxval,len = N_colLevels+1)
  } else if (discreteLevels == TRUE) {
    if (is.null(colLevels)){
      N_colLevels <- length(unique(theme[!is.nan(theme)]))
      Breakpoints <- c(sort(unique(theme[!is.nan(theme)])),2*max(theme[!is.nan(theme)]))
    } else {N_colLevels <- length(colLevels) - 1
    Breakpoints <- colLevels}}
  
  if (typeof(colPalette)=="closure") {
    colPalette <- colPalette(N_colLevels)
  } else if (typeof(colPalette)=="character") {
    if (length(colPalette) < N_colLevels){
      stop(sprintf('Length of colPalette (%d) is lower than number of colors (%d).',length(colPalette),N_colLevels))
    }
    colPalette <- colPalette[1:N_colLevels] }
  
  if (length(theme)==OCN$RN$nNodes && (length(theme)==OCN$AG$nNodes)){
    if (isTRUE(chooseAggregation == "RN")){
      byRN = TRUE
    } else if (isTRUE(chooseAggregation == "AG")){
      byRN = FALSE
    } else {
      byRN = FALSE
      #stop('RN$nNodes = AG$nNodes, and chooseAggregation has not been specified.')
    }
  } else if (length(theme)==OCN$RN$nNodes){
    byRN = TRUE
  } else if (length(theme)==OCN$AG$nNodes){
    byRN = FALSE
  } else {
    stop('theme has invalid length')
  }
  
  if (length(cex)>1 && length(cex) != length(theme)){
    stop('cex has invalid length')
  }
  
  if (length(pch)>1 && length(pch) != length(theme)){
    stop('pch has invalid length')
  }
  
  if (chooseCM==TRUE && is.logical(chooseCM)){
    chooseCM <- which(OCN$CM$A==max(OCN$CM$A))
  } else if (isFALSE(chooseCM)) {
    chooseCM <- 1:OCN$nOutlet
  }
  
  if (exactDraw==TRUE){
    X <- OCN$FD$XDraw[which(OCN$FD$toRN!=0)] # at RN level
    Y <- OCN$FD$YDraw[which(OCN$FD$toRN!=0)] # at RN level
    Xc <- OCN$CM$XContourDraw
    Yc <- OCN$CM$YContourDraw
  } else {
    X <- OCN$RN$X
    Y <- OCN$RN$Y
    Xc <- OCN$CM$XContour
    Yc <- OCN$CM$YContour
  }
  
  if (length(cex)==1){
    cex_vec <- cex*rep(1,length(theme))
  } else {cex_vec <- cex}
  
  if (length(pch)==1){
    pch_vec <- pch*rep(1,length(theme))
  } else {pch_vec <- pch}
  
  AvailableNodes <- setdiff(which(OCN$RN$toCM %in% chooseCM),OCN$RN$outlet)

  if (is.null(dots$xlim)){ dots$xlim <- c(min(X[OCN$RN$toCM %in% chooseCM]),max(X[OCN$RN$toCM %in% chooseCM])) }
  if (is.null(dots$ylim)){ dots$ylim <- c(min(Y[OCN$RN$toCM %in% chooseCM]),max(Y[OCN$RN$toCM %in% chooseCM])) }
  dots$x <- dots$xlim
  dots$y <- dots$ylim
  
  if (!add) {do.call(plot, dots)}
              
  xy_lim <- par("usr")
  
  # plot(c(min(X[OCN$FD$toCM %in% chooseCM]),max(X[OCN$FD$toCM %in% chooseCM])),
  #      c(min(Y[OCN$FD$toCM %in% chooseCM]),max(Y[OCN$FD$toCM %in% chooseCM])),
  #      type="n",xlab=" ",ylab=" ",axes=FALSE,asp=1)
  
  if (!is.null(backgroundColor)){
    if ((length(chooseCM) > 1) && (length(backgroundColor)==1) ){
      backgroundColor=rep(backgroundColor,length(chooseCM))}
    k <- 1
    for (i in chooseCM){
      for (j in 1:length(Xc[[i]])){
        polygon(Xc[[i]][[j]],Yc[[i]][[j]],col=backgroundColor[k],lty=0)
      }
      k <- k + 1
    }
  }
  
  for (i in AvailableNodes){
    #rn <- OCN$FD$toRN[i]
    reach <- OCN$RN$toAGReach[i]
    if (OCN$RN$A[i]>=OCN$thrA & 
        abs(X[i]-X[OCN$RN$downNode[i]]) <= 1.001*OCN$cellsize & 
        abs(Y[i]-Y[OCN$RN$downNode[i]]) <= 1.001*OCN$cellsize &
        X[i] >= xy_lim[1]-1*OCN$cellsize &
        X[i] <= xy_lim[2]+1*OCN$cellsize &
        Y[i] >= xy_lim[3]-1*OCN$cellsize & 
        Y[i] <= xy_lim[4]+1*OCN$cellsize) {
      if ( (byRN==TRUE && (is.nan(theme[i])==TRUE | is.na(theme[i])==TRUE)) ||
           (byRN==FALSE && (is.nan(theme[reach])==TRUE | is.na(theme[reach])==TRUE)) ||
           (byRN==TRUE && cutoff==TRUE && (theme[i] < min(Breakpoints) || theme[i] > max(Breakpoints))) ||
           (byRN==FALSE && cutoff==TRUE && (theme[reach] < min(Breakpoints) || theme[reach] > max(Breakpoints))) )  {
        hexcolor <- nanColor
      } else {
        if (byRN==TRUE){val <- theme[i]} else {val <- theme[reach]}
        
        colvalue <- which(Breakpoints > val)[1] - 1  
        if (isTRUE(colvalue==0)) {colvalue <- 1}
        if (is.na(colvalue)) {colvalue <- N_colLevels}
        
        hexcolor <- colPalette[colvalue]
      }
      if (drawNodes==TRUE | all(is.na(theme))){
        lines(c(X[i],X[OCN$RN$downNode[i]]),c(Y[i],Y[OCN$RN$downNode[i]]),
              lwd=min_lwd+(max_lwd-min_lwd)*(OCN$RN$A[i]/max(OCN$RN$A[AvailableNodes]))^0.5,col=riverColor)
      } else {
        lines(c(X[i],X[OCN$RN$downNode[i]]),c(Y[i],Y[OCN$RN$downNode[i]]),
              lwd=min_lwd+(max_lwd-min_lwd)*(OCN$RN$A[i]/max(OCN$RN$A[AvailableNodes]))^0.5,col=hexcolor)}
    }
  }
  
  if (drawNodes==TRUE && exactDraw==FALSE){
    if (byRN==TRUE){
      nodes <- which(OCN$RN$toCM %in% chooseCM)
    } else {nodes <- which(OCN$AG$toCM %in% chooseCM)}
    for (i in nodes){
      if (is.nan(theme[i]) || (cutoff==TRUE && (theme[i] < min(Breakpoints) || theme[i] > max(Breakpoints)) )){
        hexcolor <- nanColor
      } else {
        colvalue <- which(Breakpoints > theme[i])[1] - 1
        if (isTRUE(colvalue==0)) {colvalue <- 1}
        if (is.na(colvalue)) {colvalue <- N_colLevels}
        #colvalue <- 1+round((N_colLevels-1)*max(0,min(1,(theme[i]-minval)/(maxval-minval))))
        hexcolor <- colPalette[colvalue]
      }
      if (byRN==TRUE){
        points(OCN$RN$X[i],OCN$RN$Y[i],bg=hexcolor,pch=pch_vec[i],cex=cex_vec[i])
      } else if (nodeType=="upstream") {
        points(OCN$AG$X[i],OCN$AG$Y[i],bg=hexcolor,pch=pch_vec[i],cex=cex_vec[i])  
      } else if (nodeType=="downstream") {
        points(OCN$AG$XReach[i],OCN$AG$YReach[i],bg=hexcolor,pch=pch_vec[i],cex=cex_vec[i])  
      }
    }
  }
  
  if (drawNodes==TRUE && exactDraw==TRUE){
    if (byRN==TRUE){
      nodes <- which(OCN$RN$toCM %in% chooseCM)
    } else {nodes <- which(OCN$AG$toCM %in% chooseCM)}
    for (i in nodes){
      if (is.nan(theme[i]) || (cutoff==TRUE && (theme[i] < min(Breakpoints) || theme[i] > max(Breakpoints)) )){
        hexcolor <- nanColor
      } else {
        colvalue <- which(Breakpoints > theme[i])[1] - 1  
        if (isTRUE(colvalue==0)) {colvalue <- 1}
        if (is.na(colvalue)) {colvalue <- N_colLevels}
        #colvalue <- 1+round((N_colLevels-1)*max(0,min(1,(theme[i]-minval)/(maxval-minval))))
        hexcolor <- colPalette[colvalue]
      }
      if (byRN==TRUE){
        node <- 1:OCN$RN$nNodes #which(OCN$FD$X==OCN$RN$X[i] & OCN$FD$Y==OCN$RN$Y[i])
      } else if (byRN==FALSE && nodeType=="upstream") {
        node <- which(OCN$RN$X==OCN$AG$X[i] & OCN$RN$Y==OCN$AG$Y[i])
      } else if (byRN==FALSE && nodeType=="downstream") {
        node <- which(OCN$RN$X==OCN$AG$XReach[i] & OCN$RN$Y==OCN$AG$YReach[i])
      }
      points(X[node],Y[node],bg=hexcolor,pch=pch_vec[i],cex=cex_vec[i])  
    }
  }
  
  
  if (addLegend) {
    if (discreteLevels==FALSE){
      if (is.null(args_imagePlot$smallplot)){
      args_imagePlot$smallplot <- c(0.88, 0.9,par()$plt[3],par()$plt[4])}
      if (is.null(args_imagePlot$col)){args_imagePlot$col <- colPalette} 
      if (is.null(args_imagePlot$legend.only)){args_imagePlot$legend.only <- TRUE}
      if (is.null(args_imagePlot$zlim)){args_imagePlot$zlim <- c(minval, maxval)}
      do.call(imagePlot, args_imagePlot)
      # imagePlot(col=colPalette,legend.only=TRUE,zlim=c(minval,maxval),
      #            smallplot=smallplot)
    } else {
      if (is.null(colLevels)){
        str <- NULL
        for (level in 1:N_colLevels){
          str <- c(str, as.character(round(1000*Breakpoints[level])/1000) )}
      } else { 
        str <- vector(mode="character", N_colLevels)
        for (level in 1:(N_colLevels-1)){
          str[level] <- paste("[",as.character(round(1000*Breakpoints[level])/1000),"; ",
                              as.character(round(1000*Breakpoints[level+1])/1000),")",sep="")
        } 
        str[N_colLevels] <- paste("[",as.character(round(1000*Breakpoints[N_colLevels])/1000),"; ",
                                  as.character(round(1000*Breakpoints[N_colLevels+1])/1000),"]",sep="")
      }
      
      if (is.null(args_legend$x)){args_legend$x <- 1.01*max(X[OCN$RN$toCM %in% chooseCM])}
      if (is.null(args_legend$y)){args_legend$y <- max(Y[OCN$RN$toCM %in% chooseCM])}
      if (is.null(args_legend$legend)){args_legend$legend <- str}
      if (is.null(args_legend$fill)){args_legend$fill <- colPalette}
      if (is.null(args_legend$ncol)){args_legend$ncol <- ceiling(N_colLevels/20)}
      if (is.null(args_legend$xpd)){args_legend$xpd <- TRUE}
      if (is.null(args_legend$cex)){args_legend$cex <- 0.8}
      if (is.null(args_legend$bty)){args_legend$bty <- "n"}
      
      do.call(legend,args_legend)
      # legend(x=1.01*max(X[OCN$RN$toCM %in% chooseCM]),y=max(Y[OCN$RN$toCM %in% chooseCM]),
      #        str,fill=colPalette,ncol=ceiling(N_colLevels/20), xpd=TRUE, cex=0.8, bty="n")
    }
  }
  invisible()
}