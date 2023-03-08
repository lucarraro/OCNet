
draw_subcatchments_OCN <- function(OCN,
                                   theme = NULL,
                                   drawRiver = TRUE,
                                   colPalette = NULL,
                                   colLevels = NULL,
                                   riverColor = NULL,
                                   addLegend = NULL,
                                   min_lwd = 0.5,
                                   max_lwd = 5,
                                   add = FALSE,
                                   args_imagePlot = list(),
                                   ...){
  
  dots <- list(...)
  if (is.null(dots$axes)){
    if ("xlim" %in% names(dots) | "ylim" %in% names(dots)) {dots$axes <- TRUE} else {dots$axes <- FALSE}}
  if (is.null(dots$xlab)){dots$xlab <- ""}
  if (is.null(dots$ylab)){dots$ylab <- ""}
  if (is.null(dots$asp)){dots$asp <- 1}
  
  if (is.null(dev.list()) & add==TRUE){
    add <- FALSE
    warning("'add' will be ignored as there is no existing plot")
  }
  
  if (length(OCN$RN$nNodes)==0){
    stop('Missing fields in OCN. You should run aggregate_OCN prior to draw_subcatchments_OCN.')
  }
  
  if (is.null(colPalette)){
    if (is.null(theme)){
      colPalette <- c("#009900", # green
                               "#FFFF00", # yellow
                               "#FF9900", # orange
                               "#FF0000", # red
                               "#FF00FF", # fuchsia
                               "#9900CC", # violet
                               "#555555", # grey 30%
                               "#BBBBBB") # grey 70%
    } else { colPalette=colorRampPalette(c("yellow","red","black")) }
  } else if (typeof(colPalette)=="closure" & is.null(theme)) { colPalette <- colPalette(8)}
  
  if (!is.null(theme)){
    if (length(colLevels)<3){N_colLevels <- 1000} else {N_colLevels <- colLevels[3]}
    if (typeof(colPalette)=="closure"){
      colPalette <- colPalette(N_colLevels)
    } else {
      if (length(colPalette) < N_colLevels){ stop(sprintf('Length of colPalette (%d) is lower than number of colors (%d).',length(colPalette),N_colLevels)) }
      colPalette <- colPalette[1:N_colLevels]}
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
  }
  
  if (is.null(riverColor)){
    if (is.null(theme)){riverColor <- "black"}
    else riverColor <- "#00BFFF"
  }
  
  if (is.null(theme)){
    ## Greedy algorithm for coloring subcatchment map
    ## Create list of nodes for the greedy algorithm
    W <- OCN$SC$W
    ListNodes <- numeric(OCN$SC$nNodes)
    for (i in 1:(OCN$SC$nNodes)){
      degree <- colSums(W)
      degree <- degree[degree>0]
      if (length(degree)>0){
        MinDegree <- which(colSums(W)==min(degree))
        NodeRemoved <- MinDegree[1]
        ListNodes[OCN$SC$nNodes+1-i] <- NodeRemoved
        W[,NodeRemoved] <- 0
        W[NodeRemoved,] <- 0
      } else {
        MissingNode <- setdiff(1:OCN$SC$nNodes,ListNodes)[1]
        ListNodes[OCN$SC$nNodes+1-i] <- MissingNode
      }
    }
    
    ## Attribute colors
    ColorList <- 1
    ColorID <- numeric(OCN$SC$nNodes)
    for (i in 1:OCN$SC$nNodes){
      connected_nodes <- which(OCN$SC$W[ListNodes[i],1:OCN$SC$nNodes]==1)
      k <- setdiff(ColorList,ColorID[connected_nodes])[1]
      if (is.na(k)){
        ColorList <- c(ColorList,max(ColorList)+1)
        k <- max(ColorList)
      }
      ColorID[ListNodes[i]] <- k
    }
    colPalette <- colPalette[ColorList]
  }
  
  ## plot subcatchment map
  kol <- numeric(OCN$FD$nNodes)
  for (k in 1:OCN$SC$nNodes){
    if (is.null(theme)) {kol[OCN$SC$toFD[[k]]] <- ColorID[k]}
    else { colvalue <- which(Breakpoints > theme[k])[1] - 1  
    if (isTRUE(colvalue==0)) {colvalue <- 1}
    if (is.na(colvalue)) {colvalue <- N_colLevels}
    kol[OCN$SC$toFD[[k]]] <- colvalue
    }
  }
  
  if (OCN$FD$nNodes < OCN$dimX*OCN$dimY){
    if (isTRUE(OCN$typeInitialState=="custom")){
      Color_SC <- matrix(NaN,OCN$dimY,OCN$dimX)
      Color_SC[OCN$FD$toDEM] <- kol
      Color_SC <- Color_SC[seq(OCN$dimY,1,-1), ]
    } else { # real river
      Color_SC <- matrix(NaN,OCN$dimX,OCN$dimY)
      Color_SC[OCN$FD$toDEM] <- kol
      Color_SC <- Color_SC[,seq(OCN$dimY,1,-1)]
      Color_SC <- t(Color_SC)
    }
  } else {Color_SC <- matrix(data=kol,nrow=OCN$dimY,ncol=OCN$dimX)}
  
  
  if (length(OCN$xllcorner)==0){xllcorner <- min(OCN$FD$X)[1]} else {xllcorner <- OCN$xllcorner}
  if (length(OCN$yllcorner)==0){yllcorner <- min(OCN$FD$Y)[1]} else {yllcorner <- OCN$yllcorner}
  
  Xvec <- seq(xllcorner,xllcorner+(OCN$dimX-1)*OCN$cellsize,OCN$cellsize)
  Yvec <- seq(yllcorner,yllcorner+(OCN$dimY-1)*OCN$cellsize,OCN$cellsize)
  
  #old.par <- par(no.readonly = TRUE)
  #on.exit(par(old.par))
  #par(bty="n")
  dots$x <- Xvec; dots$y <- Yvec; dots$z <- t(Color_SC); dots$col <- colPalette
  dots$xlim <- range(OCN$FD$X); dots$ylim <- range(OCN$FD$Y)
  dots$add <- add
  
  do.call(image,dots)
  
  # image(Xvec, Yvec,
  #       t(Color_SC),col=colPalette,xlab=" ",ylab=" ",asp=1, axes=FALSE) #
  # attributing colors in reverse order should increase overall contrast
  
  if (drawRiver==TRUE){
    ## plot OCN
    AvailableNodes <- setdiff(1:OCN$RN$nNodes,OCN$RN$outlet)
    #points(OCN$FD$X[OCN$FD$outlet],OCN$FD$Y[OCN$FD$outlet],pch=22,col="#000000",bg="#000000")
    
    for (i in AvailableNodes){
      if (OCN$RN$A[i]>=OCN$thrA  & abs(OCN$RN$X[i]-OCN$RN$X[OCN$RN$downNode[i]])<=1.001*OCN$cellsize & 
          abs(OCN$RN$Y[i]-OCN$RN$Y[OCN$RN$downNode[i]])<=1.001*OCN$cellsize  ) {
        lines(c(OCN$RN$X[i],OCN$RN$X[OCN$RN$downNode[i]]),c(OCN$RN$Y[i],OCN$RN$Y[OCN$RN$downNode[i]]),
              lwd=min_lwd + (max_lwd-min_lwd)*(OCN$RN$A[i]/(OCN$FD$nNodes*OCN$cellsize^2))^0.5,
              col=riverColor)}
    }}
  
  if (is.null(addLegend) & !is.null(theme)){addLegend <- TRUE} else {addLegend <- FALSE}
  if (addLegend) {
    if (is.null(args_imagePlot$smallplot)){
      args_imagePlot$smallplot <- c(0.88, 0.9,par()$plt[3],par()$plt[4])}
    if (is.null(args_imagePlot$col)){args_imagePlot$col <- colPalette} 
    if (is.null(args_imagePlot$legend.only)){args_imagePlot$legend.only <- TRUE}
    if (is.null(args_imagePlot$zlim)){args_imagePlot$zlim <- c(minval, maxval)}
    do.call(imagePlot, args_imagePlot)
    # imagePlot(col=colPalette,legend.only=TRUE,zlim=c(minval,maxval),
    #           smallplot=c(0.88, 0.9,par()$plt[3],par()$plt[4]))
  }
  invisible()
}