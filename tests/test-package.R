#
# Package tests
#
#
# 

require(OCNet)

# agg contains length, and cumsum of all numeric entries
ocn_test <- function(obj, agg, tol=1e-6, label="") {
  
  stopifnot( length(obj)==agg[1])
  
  tmp <- sum( unlist( obj[unlist( lapply(obj, is.numeric))]))
  stopifnot( abs(tmp-agg[2]) < tol)
  
  tmp <- sum( unlist( obj$FD[unlist( lapply(obj$FD, is.numeric))]))
  stopifnot( abs(tmp-agg[3]) < tol)

  if (!is.null(obj$CR)){
    tmp <- sum( unlist( obj$CR[unlist( lapply(obj$CR, is.numeric))]))
    stopifnot( abs(tmp-agg[4]) < tol)
  }
  return(TRUE)    
}

get_agg <- function( obj) {
  agg <-
    c( length(obj),
       sum( unlist( obj[unlist( lapply(obj, is.numeric))])),
       sum( unlist( obj$FD[unlist( lapply(obj$FD, is.numeric))]))
    )
  if (!is.null(obj$CR)){
    agg[4] <-  sum( unlist( obj$CR[unlist( lapply(obj$CR, is.numeric))]))
  }
  print( agg, digits=16)
  return(agg)
}

ocn <- create_OCN(12,3,  seed=1)
ocn_test( ocn, c(12.0, 1458.5, 1095.))
ocn <- create_OCN(12, 3, seed=1, saveEnergy = TRUE)
ocn_test( ocn, c( 13,93892.78110789257,  1095))
  
  
# the following should cause errors!
try( create_OCN(12, -3))
try( create_OCN(4, 4, outletSide = "K"))


# To be completed....
