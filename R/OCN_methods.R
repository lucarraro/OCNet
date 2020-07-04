# methods handling for the package.

# Naming: 
# OCN
# OCNlandscape
# OCNaggregate
# print, summary and plot
#
#  Classes with inheritance:			
#    https://www.datamentor.io/r-programming/inheritance/
#  OCN > OCNlandscape > OCNaggregate
#


print.OCN <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
  cat(paste0("Optimal channel network of dimension ",x$dimX,"x",x$dimX,".\n"))
    
  invisible(x)
}

summary.OCN <- function(object, ...) {
  cat("Object of class `OCN`;\n   to be feed to function `landscape_OCN` for further processing.")
  cat(paste0("\nOCNet of dimension: ",object$dimX,"x",object$dimX," with ",object$nOutlet," outlets."))
  cat(paste0("\nCooling scheme ... \n"))
  invisible(object)
}

plot.OCN <- function(x, ...) {        
    draw_simple_OCN( x, ...)
    invisible()  
}


##########################################################
print.OCNlandscape <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call_landscape), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
  cat(paste0("Landscaped optimal channel network of dimension ",x$dimX,"x",x$dimX,".\n"))    
  invisible(x)
}

summary.OCNlandscape <- function(object, ...) {
  cat("Object of class `OCNlandscape`;\n   to be feed to function `aggreate_OCN` for further processing.")
  cat(paste0("\nOCNet of dimension: ",object$dimX,"x",object$dimX," with ",object$nOutlet," outlets."))
  cat(paste0("\nCooling scheme ... \n"))
  invisible(object)
}

plot.OCNlandscape <- function(x, which=c("contour","elev2D","elev3D","elev3Drgl","simple"), ...) {
  # implement which with numbers??
    type <- match.arg(which)
    switch(type,
           contour=draw_contour_OCN(x, ...),
           elev2D=draw_elev2D_OCN(x, ...),
           elev3D=draw_elev3D_OCN(x,...),
           elev3Drgl=draw_elev3Drgl_OCN(x, ...),
           simple=draw_simple_OCN(x, ...))
    invisible()  
}

##########################################################

print.OCNaggregate <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call_aggregate), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
  cat(paste0("Aggregated optimal channel network of dimension ",x$dimX,"x",x$dimX,".\n"))
    
  invisible(x)
}

summary.OCNaggregate <- function(object, ...) {
  cat("Object of class `OCN`;\n   to be feed to function `landscape_OCN` for further processing.")
  cat(paste0("\nOCNet of dimension: ",object$dimX,"x",object$dimX," with ",object$nOutlet," outlets."))
  cat(paste0("\nCooling scheme ... \n"))
  invisible(object)
}

plot.OCNaggregate <- function(x, ...) {        
    draw_simple_OCN( x, ...)
    invisible()  
}

