OCN_to_igraph <- function (OCN, 
                           level) {
  
  if (missing(OCN)) {
    stop("Input OCN cannot be missing")
  }
  if (!(level %in% c("FD","RN","AG"))) {
    stop("Invalid level")
  }
  if (!("RN" %in% names(OCN)) && (level %in% c("RN","AG"))){
    stop('Missing aggregation level in OCN. Run landscape_OCN and/or aggregate_OCN prior to OCN_to_SSN.')
  }
  
  sub_OCN <- NULL
  eval(parse(text=(paste("sub_OCN <- OCN$",level,sep=""))))
  
  mm <- as.dgCMatrix.spam(sub_OCN$W)
  g <- graph_from_adjacency_matrix(mm)
  
  return(g)
}
  