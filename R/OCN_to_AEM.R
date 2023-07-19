OCN_to_AEM <- function(OCN, level="AG", weight = NULL,
                       resistance = "length", moranI = FALSE) {

  if (length(OCN$RN$nNodes)==0){
    stop('Missing aggregation level in OCN. Run landscape_OCN and/or aggregate_OCN prior to OCN_to_AEM.')
  }
  
  if (resistance=="time" & length(OCN[[level]]$velocity)==0){
    stop('Missing velocities. resistance = "time" cannot be used.
         Use OCNet::rivergeometry_OCN or rivnet::hydro_river to compute velocities.')
  }

  se.mat.OCN <- matrix(0,OCN[[level]]$nNodes, OCN[[level]]$nNodes)
  for (i in 1:OCN[[level]]$nNodes) se.mat.OCN[OCN[[level]]$upstream[[i]],i] <- 1

  edges.OCN <- cbind(OCN[[level]]$downNode, 1:OCN[[level]]$nNodes)
  colnames(edges.OCN) <- c("from","to")
  bin.mat.OCN <- list(se.mat=se.mat.OCN, edges=edges.OCN)

  if (resistance=="length"){
    resistance <- OCN[[level]]$leng
  } else if (resistance=="time"){
    resistance <- OCN[[level]]$leng/OCN[[level]]$velocity
  } else {stop("Invalid resistance")}

  if (is.null(weight)){
    weight <- 1+numeric(OCN[[level]]$nNodes)
  } else if (typeof(weight)=="closure"){
    weight <- weight(resistance)
  } else if (weight=="gravity"){
    weight <- max(resistance)/resistance
  } else if (weight=="linear"){
    weight <- 1 - (resistance/max(resistance))
  } else if (weight=="exponential"){
    weight <- exp(-resistance/max(resistance))
  } else if (weight=="parabolic"){
    weight <- 1 - (resistance/max(resistance))^2
  }  else {stop("Invalid weight")}
  weight[weight==Inf] <- 0 # it doesn't matter which value

  res <- aem(bin.mat.OCN, weight=weight)
  
  if (moranI){
    matW <- t(as.matrix.spam(OCN[[level]]$W * weight))
    listW <- suppressWarnings(mat2listw(matW))
    moranI <- moran.randtest(res$vectors, listW)
    res[["moranI"]] <- moranI
  }
  
  invisible(res)
}
