
rivergeometry_OCN <- function(OCN,
                              widthMax=1,
                              depthMax=1,
                              velocityMax=1,
                              expWidth=NaN,
                              expDepth=NaN,
                              expVelocity=NaN){
  
  exponents=c(expWidth,expDepth,expVelocity)

  if (any(exponents<0,na.rm=TRUE) || any(exponents>1,na.rm=TRUE)){
    stop('Invalid exponent. expWidth, expDepth, expVelocity must be bound between 0 and 1.')
  }

  if (sum(is.nan(exponents))==3){
    expWidth <- 0.5
    expDepth <- 0.4
    expVelocity <- 0.1
  } else if (sum(is.nan(exponents))==2) {
    stop('You must specify at least two inputs among expWidth, expDepth, expVelocity.')
  } else if (sum(is.nan(exponents))==1) {
    if (sum(exponents[!is.nan(exponents)])>1) {
      stop('expWidth + expDepth + expVelocity must be equal to 1.')
    } else {exponents[is.nan(exponents)] = 1 - sum(exponents[!is.nan(exponents)])}
  } else if (sum(exponents) != 1){
    stop('expWidth + expDepth + expVelocity must be equal to 1.')
  }

  if (!("AG" %in% names(OCN))){
    stop('Missing fields in OCN. You should run aggregate_OCN prior to rivergeometry_OCN.')
  }

  width_RN <- widthMax*(OCN$RN$A/(OCN$FD$nNodes*OCN$cellsize^2))^expWidth
  depth_RN <- depthMax*(OCN$RN$A/(OCN$FD$nNodes*OCN$cellsize^2))^expDepth
  velocity_RN <- velocityMax*(OCN$RN$A/(OCN$FD$nNodes*OCN$cellsize^2))^expVelocity
  
  A_AG <- 0.5*(OCN$AG$A + OCN$AG$AReach)
  
  width_AG <- widthMax*(A_AG/(OCN$FD$nNodes*OCN$cellsize^2))^expWidth
  depth_AG <- depthMax*(A_AG/(OCN$FD$nNodes*OCN$cellsize^2))^expDepth
  velocity_AG <- velocityMax*(A_AG/(OCN$FD$nNodes*OCN$cellsize^2))^expVelocity

  # copy variables into OCN list
  OCN$RN[["width"]] <- width_RN
  OCN$RN[["depth"]] <- depth_RN
  OCN$RN[["velocity"]] <- velocity_RN
  OCN$AG[["width"]] <- width_AG
  OCN$AG[["depth"]] <- depth_AG
  OCN$AG[["velocity"]] <- velocity_AG

  OCN$widthMax <- widthMax
  OCN$depthMax <- depthMax
  OCN$velocityMax <- velocityMax
  
  OCN$expWidth <- expWidth
  OCN$expDepth <- expDepth
  OCN$expVelocity  <- expVelocity 

  invisible(OCN)
}