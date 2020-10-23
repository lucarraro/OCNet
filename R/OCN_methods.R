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
  
  writeLines(strwrap( paste0("Optimal channel network of dimension ",x$dimX,"x",x$dimX,
    {if (x$nOutlet==1) " with one single outlet" else
      if (x$nOutlet<(2*(x$dimX + x$dimY - 2))) paste0(" with ",x$nOutlet," outlets")},
    {if (x$periodicBoundaries) " and periodic boundaries"}, ".")
  ))
  cat("\n")

  invisible(x)
}

summary.OCN <- function(object, ...) {
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  writeLines(strwrap( paste0("Optimal channel network of dimension ",object$dimX,"x",object$dimX,
        {if (object$nOutlet==1) " with one single outlet" else
         if (object$nOutlet<(2*(object$dimX + object$dimY - 2))) paste0(" with ",object$nOutlet," outlets")},
       {if (object$periodicBoundaries) " and periodic boundaries"}, ". Size of a pixel is ",object$cellsize," planar units.")))
  
  writeLines(strwrap( paste0("Cooling parameters were 'coolingRate'=", 
    object$coolingRate," and 'initialNoCoolingPhase'=",object$initialNoCoolingPhase,
                             " with exponent ",object$expEnergy," for ",object$nIter," iterations.")
  ))
  
  cat("\nObject of class `OCN`;\n   to be feed to function `landscape_OCN` for further processing.\n\n")
                      
  
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
  
  writeLines(strwrap( paste0("Landscaped optimal channel network of dimension ",x$dimX,"x",x$dimX,
                             {if (x$nOutlet==1) " with one single outlet" else
                               if (x$nOutlet<(2*(x$dimX + x$dimY - 2))) paste0(" with ",x$nOutlet," outlets")},
                             {if (x$periodicBoundaries) " and periodic boundaries"}, ".")
  ))
  cat("\n")
  
  invisible(x)
}

summary.OCNlandscape <- function(object, ...) {
  
  cat("\nCall:\n", paste(deparse(object$call_landscape), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  writeLines(strwrap( paste0("Optimal channel network of dimension ",object$dimX,"x",object$dimX,
                             {if (object$nOutlet==1) " with one single outlet" else
                               if (object$nOutlet<(2*(object$dimX + object$dimY - 2))) paste0(" with ",object$nOutlet," outlets")},
                             {if (object$periodicBoundaries) " and periodic boundaries"}, ". Size of a pixel is ",object$cellsize," planar units.")))
  
  writeLines(strwrap( paste0("Cooling parameters were 'coolingRate'=", 
                             object$coolingRate," and 'initialNoCoolingPhase'=",object$initialNoCoolingPhase,
                             " with exponent ",object$expEnergy," for ",object$nIter," iterations.")
  ))
  # periodic boundary, total area in units,  elevation range, slope of outlet.
  #  elevation adjustment of multiple outlets, arg=optimizeDZ,
  
  writeLines(strwrap( paste0("The landscape covers ",
            object$dimX*object$dimY*object$cellsize," planar units and stretches over ", 
            round(diff( range(object$FD$slope)),2), " elevation units (range: ",
            round(min(object$FD$slope),2),"-",round( max(object$FD$slope),2),"). Slope of outlet is ",
            object$slope0,". ",
            if(object$nOutlet>1) {
              ifelse(object$optimizeDZ, 
                     "Differences at the catchment borders have beend minimized.",
                     "No elevation adjustments at the catchment borders have been performed.")} )
                      ))
  
  
  cat("\nObject of class `OCNlandscape`;\n   to be feed to function `aggreate_OCN` for further processing.\n\n")
  
  invisible(object)
}



plot.OCNlandscape <- function(x, which=c("contour","simple","elev2D","elev3D","elev3Drgl"), ...) {

  if (is.numeric(which)){
    if (any(which < 1) || any(which > 5)) 
      stop("'which' must be in 1:5 or in c('contour','simple','elev2D','elev3D','elev3Drgl')")
    which <- c("contour","simple","elev2D","elev3D","elev3Drgl")[which]
  } # else which <- tolower(which)
  # the to lower does not work well. Want to keep the capital 'D'
  type <- match.arg(which)
  switch(type,
         contour=draw_contour_OCN(x, ...),
         elev2D=draw_elev2D_OCN(x, ...),
         elev3D=draw_elev3D_OCN(x, ...),
         elev3Drgl=draw_elev3Drgl_OCN(x, ...),
         simple=draw_simple_OCN(x, ...))
  invisible()  
}

##########################################################

print.OCNaggregate <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call_aggregate), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")

  writeLines(strwrap( paste0("Aggregated optimal channel network of dimension ",x$dimX,"x",x$dimX,
                             {if (x$nOutlet==1) " with one single outlet" else
                               if (x$nOutlet<(2*(x$dimX + x$dimY - 2))) paste0(" with ",x$nOutlet," outlets")},
                             {if (x$periodicBoundaries) " and periodic boundaries"}, ".")
  ))
  cat("\n")
  
  invisible(x)
}

summary.OCNaggregate <- function(object, ...) {
   # all of the above, number of nodes at RNlevel, AGlevel, at subcatchment level
  cat("\nCall:\n", paste(deparse(object$call_aggregate), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  writeLines(strwrap( paste0("Optimal channel network of dimension ",object$dimX,"x",object$dimX,
                             {if (object$nOutlet==1) " with one single outlet" else
                               if (object$nOutlet<(2*(object$dimX + object$dimY - 2))) paste0(" with ",object$nOutlet," outlets")},
                             {if (object$periodicBoundaries) " and periodic boundaries"}, ". Size of a pixel is ",object$cellsize," planar units.")))
  
  writeLines(strwrap( paste0("Cooling parameters were 'coolingRate'=", 
                             object$coolingRate," and 'initialNoCoolingPhase'=",object$initialNoCoolingPhase,
                             " with exponent ",object$expEnergy," for ",object$nIter," iterations.")
  ))

  writeLines(strwrap( paste0("The landscape covers ",
                             object$dimX*object$dimY*object$cellsize," planar units and stretches over ", 
                             round(diff( range(object$FD$slope)),2), " elevation units (range: ",
                             round(min(object$FD$slope),2),"-",round( max(object$FD$slope),2),"). Slope of outlet is ",
                             object$slope0,". ",
                             if(object$nOutlet>1) {
                               ifelse(object$optimizeDZ, 
                                      "Differences at the catchment borders have beend minimized.",
                                      "No elevation adjustments at the catchment borders have been performed.")} )
  ))

  cat('\n')
  
  writeLines(strwrap( paste0("The aggregated landscape covers ", 
    object$RN$nNodes, " nodes at river network, ",
    object$AG$nNodes, " nodes at aggregated network and ",
    object$SC$nNodes, " nodes at subcatchment level.")))
  cat('\n')
  
  invisible(object)
}

plot.OCNaggregate <- function(x, which=c("thematic","subcatchment","contour","simple",
                                         "elev2D","elev3D","elev3Drgl"), 
                              theme=NULL, ...) {     
  if (is.numeric(which)){
    if (any(which < 1) || any(which > 7)) 
      stop("'which' must be in 1:7 or in c('thematic','subcatchment','contour','simple','elev2D','elev3D','elev3Drgl')")
    which <- c("thematic","subcatchment","contour","simple","elev2D","elev3D","elev3Drgl")[which]
  } # else which <- tolower(which)
  # the to lower does not work well. Want to keep the capital 'D'
  type <- match.arg(which)
  
  if (type=='thematic' && is.null(theme)) theme <- x$AG$A  #  might need some adjustment!!!!!

  switch(type,
         thematic=draw_thematic_OCN(theme, x, ...),  
            #  might need some adjustment!!!!!
         subcatchment=draw_subcatchments_OCN(x, ...),
         contour=draw_contour_OCN(x, ...),
         elev2D=draw_elev2D_OCN(x, ...),
         elev3D=draw_elev3D_OCN(x, ...),
         elev3Drgl=draw_elev3Drgl_OCN(x, ...),
         simple=draw_simple_OCN(x, ...))
  
  invisible()  
}
