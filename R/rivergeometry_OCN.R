
rivergeometry_OCN <- function(OCN,
                              width_max=1,
                              depth_max=1,
                              velocity_max=1,
                              exp_width=NaN,
                              exp_depth=NaN,
                              exp_velocity=NaN){
  
  exponents=c(exp_width,exp_depth,exp_velocity)

  if (any(exponents<0,na.rm=TRUE) || any(exponents>1,na.rm=TRUE)){
    stop('Invalid exponent. exp_width, exp_depth, exp_velocity must be bound between 0 and 1.')
  }

  if (sum(is.nan(exponents))==3){
    exp_width <- 0.5
    exp_depth <- 0.4
    exp_velocity <- 0.1
  } else if (sum(is.nan(exponents))==2) {
    stop('You must specify at least two inputs among exp_width, exp_depth, exp_velocity.')
  } else if (sum(is.nan(exponents))==1) {
    if (sum(exponents[!is.nan(exponents)])>1) {
      stop('exp_width + exp_depth + exp_velocity must be equal to 1.')
    } else {exponents[is.nan(exponents)] = 1 - sum(exponents[!is.nan(exponents)])}
  } else if (sum(exponents) != 1){
    stop('exp_width + exp_depth + exp_velocity must be equal to 1.')
  }

  if (!("AG" %in% names(OCN))){
    stop('Missing fields in OCN. You should run aggregate_OCN prior to rivergeometry_OCN.')
  }

  width_RN <- width_max*(OCN$RN$A/(OCN$FD$Nnodes*OCN$cellsize^2))^exp_width
  depth_RN <- depth_max*(OCN$RN$A/(OCN$FD$Nnodes*OCN$cellsize^2))^exp_depth
  velocity_RN <- velocity_max*(OCN$RN$A/(OCN$FD$Nnodes*OCN$cellsize^2))^exp_velocity
  
  width_AG <- width_max*(OCN$AG$A/(OCN$FD$Nnodes*OCN$cellsize^2))^exp_width
  depth_AG <- depth_max*(OCN$AG$A/(OCN$FD$Nnodes*OCN$cellsize^2))^exp_depth
  velocity_AG <- velocity_max*(OCN$AG$A/(OCN$FD$Nnodes*OCN$cellsize^2))^exp_velocity

  # copy variables into OCN list
  OCN$RN[["width"]] <- width_RN
  OCN$RN[["depth"]] <- depth_RN
  OCN$RN[["velocity"]] <- velocity_RN
  OCN$AG[["width"]] <- width_AG
  OCN$AG[["depth"]] <- depth_AG
  OCN$AG[["velocity"]] <- velocity_AG

  OCN$width_max <- width_max
  OCN$depth_max <- depth_max
  OCN$velocity_max <- velocity_max
  
  OCN$exp_width <- exp_width
  OCN$exp_depth <- exp_depth
  OCN$exp_velocity  <- exp_velocity 

  return(OCN)
}