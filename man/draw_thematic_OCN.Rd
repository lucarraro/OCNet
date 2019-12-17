\name{draw_thematic_OCN}
\alias{draw_thematic_OCN}
\title{Draw thematic map on an Optimal Channel Network}
\usage{
draw_thematic_OCN(theme, OCN, discreteLevels = FALSE, 
  colLevels = NULL, cutoff = FALSE,
  colPalette = colorRampPalette(c("yellow","red","black")),
  exactDraw = FALSE, chooseCM = FALSE, drawNodes = FALSE,  
  cex = 2, pch = 21, nanColor = "#0099FF",
  riverColor = "#0099FF", backgroundColor = "#999999", 
  addLegend = TRUE)
}
\arguments{
\item{theme}{Vector (of length \code{OCN$AG$Nnodes} or \code{OCN$RN$Nnodes}) expressing the spatial field of interest. 
The vector can contain \code{NA} and \code{NaN} values to identify RN or AG nodes where the theme is not defined.}

\item{OCN}{List as produced by \code{\link{aggregate_OCN}}.}

\item{discreteLevels}{Logical. If \code{FALSE}, a continuous color scheme is used. If \code{TRUE}, discrete color levels are applied. See also \code{colLevels} and examples.}

\item{colLevels}{Number of colors in the palette. If \code{discreteLevels == FALSE}, \code{colLevels} must be a vector of the form \code{c(minval, maxval, N_levels)}. The vector of breakpoints used to attribute \code{theme} values to a given color is then defined as \code{seq(minval, maxval, N_levels)}. Default is \code{minval = min(theme[!(is.nan(theme))])}, \code{maxval = max(theme[!(is.nan(theme))])}, \code{N_levels = 1000}. If \code{discreteLevels == TRUE} and \code{is.null(colLevels) == TRUE}, each unique value of \code{theme} is attributed a different color. If \code{discreteLevels == TRUE} and \code{colLevels} is a vector, \code{colLevels} is used as vector of breakpoints. In this case, the number of discrete colors is equal to \code{length(colLevels) - 1}.} 

\item{cutoff}{Logical. If \code{FALSE}, nodes whose \code{theme} value is beyond the range established by the vector of breakpoints are attributed the color corresponding to the lowest (or highest) value in the color scheme. If \code{TRUE}, such nodes are attributed the color \code{NaNcolor}.}

\item{colPalette}{Color palette used to display theme values. \code{colPalette} accepts both functions creating color palettes and vectors of colors. In the latter case, \code{length(colPalette)} must be greater than the number of color levels. See examples below and \code{\link{hcl.colors}}.}

\item{chooseCM}{Index of catchment to display (only effective if \code{OCN$N_outlet > 1}). It can be a logical or a numeric vector. If \code{FALSE}, all catchments are displayed. If \code{TRUE}, the catchment with largest area is displayed. If \code{chooseCM} is a subset of vector \code{1:length(OCN$N_outlet)}, only the catchment(s) identified by the indices in \code{chooseCM} are displayed.}

\item{exactDraw}{If \code{TRUE}, the real shape of OCNs is plotted. If flow crosses a boundary, the pixel that is not contiguous to its outlet is flipped.}

\item{drawNodes}{If \code{FALSE}, the theme is directly displayed on the river network. In this case, the edge departing from a given node is displayed with the color attributed to the node. If \code{TRUE}, the theme is displayed via markers at the locations of the nodes at the RN or AG level (depending on the length of \code{theme}).
In this case, \code{nanColor} can be used to define the color of the river network.}

\item{cex}{Only effective if \code{drawNodes == TRUE}. It sets the dimension of the markers (equivalent to 
parameter \code{cex} of function \code{\link{points}}). It can be a scalar or a vector of length \code{length(theme)}.}

\item{pch}{Only effective if \code{drawNodes == TRUE}. It sets the type of the markers (equivalent to 
parameter \code{pch} of function \code{\link{points}}). It can be a scalar or a vector of length \code{length(theme)}.}

\item{nanColor}{Color attributed to RN or AG nodes whose theme value is \code{NA} or \code{NaN}.}

\item{riverColor}{Only effective if \code{drawNodes == TRUE}. Color used to display the OCN below the nodes.}

\item{backgroundColor}{Color used in the background of the figure. It can be either a single value, or a vector with number of components 
equal to \code{length(chooseCM)}. If \code{length(backgroundColor) == length(chooseCM)}, each color is used to identify a different catchment selected in \code{chooseCM}
(corresponding to the respective outlet). If instead \code{length(chooseCM) > 1} and \code{length(backgroundColor) == 1}, all catchments are colored with the same 
\code{backgroundColor}.}

\item{addLegend}{If \code{TRUE}, add legend to the plot. If \code{discreteLevels = FALSE}, this appears a colorbar.}
}
\value{
No output is returned.
}
\description{
Function that draws OCNs with color of RN or AG nodes depending on an arbitrary theme.
}
\details{
This function can be used to show how a certain spatial field varies along the river network.
}
\examples{
# 1a) Six different ways to display contributing area at the AG level
OCN <- aggregate_OCN(landscape_OCN(OCN_20), thrA = 4)
old.par <- par(no.readonly = TRUE)
par(mfrow=c(2,3), oma = c(0, 0, 3, 0))
draw_thematic_OCN(OCN$AG$A, OCN, colPalette = hcl.colors) 
title("Continuous levels \n Colors on edges")
draw_thematic_OCN(OCN$AG$A, OCN, discreteLevels = TRUE,
		colPalette = hcl.colors) 
title("Discrete, unique levels \n Colors on edges")
draw_thematic_OCN(OCN$AG$A, OCN, discreteLevels = TRUE, 
		colLevels = c(1, 10, 50, 100, 500), 
		colPalette = hcl.colors) 
title("Discrete, user-defined levels \n Colors on edges")
draw_thematic_OCN(OCN$AG$A, OCN, drawNodes = TRUE, 
		colPalette = hcl.colors)
title("Continuous levels \n Colors on edges")
draw_thematic_OCN(OCN$AG$A, OCN,  discreteLevels = TRUE, 
		drawNodes = TRUE, colPalette = hcl.colors)
title("Discrete, unique levels \n Colors on nodes")
draw_thematic_OCN(OCN$AG$A, OCN, discreteLevels = TRUE, 
		drawNodes = TRUE, colLevels = c(1, 10, 50, 100, 500), 
		colPalette = hcl.colors)
title("Discrete, user-defined levels \n Colors on nodes")
mtext("Six different ways to display contributing area [no. pixels]", outer = TRUE, cex = 1.5)
par(old.par)

# 1b) Same as above, but use different colLevels, cutoff combinations
# with DiscreteLevels = FALSE
old.par <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
draw_thematic_OCN(OCN$AG$A, OCN, drawNodes = TRUE,
		 colLevels = c(0, 200, 1000), colPalette = hcl.colors)
title("All nodes with A > 200 pixels \n are displayed in yellow")
draw_thematic_OCN(OCN$AG$A, OCN, drawNodes = TRUE, 
		nanColor = "#00000000", colLevels = c(0, 200, 1000),
		cutoff = TRUE, colPalette = hcl.colors)
title("All nodes with A > 200 pixels \n are treated as NaN")		
par(old.par)

\donttest{
# 2) Display distance to outlet (at the RN level) along the main stem
# of an OCN
OCN <- aggregate_OCN(landscape_OCN(OCN_250_T)) # this takes some seconds
OCN <- paths_OCN(OCN, pathsRN = TRUE) # this takes some seconds

distanceToOutlet <- OCN$RN$downstreamPathLength[,OCN$RN$outlet]
farthestNode <- which(distanceToOutlet == max(distanceToOutlet))
mainStem <- OCN$RN$downstreamPath[[farthestNode]][[OCN$RN$outlet]]
theme <- rep(NaN, OCN$RN$nNodes)
theme[mainStem] <- distanceToOutlet[mainStem]

draw_thematic_OCN(theme, OCN)
title("Distance to outlet along the main stem [pixel units]")
}
}
