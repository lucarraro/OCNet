\name{draw_thematic_OCN}
\alias{draw_thematic_OCN}
\title{Draw thematic map on an Optimal Channel Network}
\usage{
draw_thematic_OCN(OCN,theme=NA*numeric(OCN$AG$nNodes),
  chooseAggregation = NULL,
  discreteLevels = FALSE, 
  colLevels = NULL, cutoff = FALSE,
  colPalette = colorRampPalette(c("yellow","red","black")),
  exactDraw = FALSE, chooseCM = FALSE, drawNodes = FALSE,  
  nodeType = "upstream",  nanColor = "#00BFFF",
  riverColor = "#00BFFF", backgroundColor = "#999999", 
  addLegend = TRUE, min_lwd = 0.5, max_lwd = 5, ...)
}
\arguments{
\item{OCN}{List as produced by \code{\link{aggregate_OCN}}. The order of arguments between \code{OCN} and \code{theme} can be swapped freely.}

\item{theme}{Vector (of length \code{OCN$AG$Nnodes} or \code{OCN$RN$Nnodes}) expressing the spatial field of interest. 
The vector can contain \code{NA} and \code{NaN} values to identify RN or AG nodes where the theme is not defined. The order of arguments between \code{OCN} and \code{theme} can be swapped freely. }

\item{chooseAggregation}{Only effective if \code{OCN$RN$nNodes == OCN$AG$nNodes}. In such case, it must be equal to either \code{"RN"} or \code{"AG"}; as a result, \code{theme} will be interpreted as a spatial field in the corresponding aggregation level.}

\item{discreteLevels}{Logical. If \code{FALSE}, a continuous color scheme is used. If \code{TRUE}, discrete color levels are applied. See also \code{colLevels} and examples.}

\item{colLevels}{Number of colors in the palette. If \code{discreteLevels == FALSE}, \code{colLevels} must be a vector of the form \code{c(minval, maxval)} or \code{c(minval, maxval, N_levels)}. The vector of breakpoints used to attribute \code{theme} values to a given color is then defined as \code{seq(minval, maxval, N_levels)}. Default is \code{minval = min(theme[!(is.nan(theme))])}, \code{maxval = max(theme[!(is.nan(theme))])}, \code{N_levels = 1000}. If \code{discreteLevels == TRUE} and \code{is.null(colLevels) == TRUE}, each unique value of \code{theme} is attributed a different color. If \code{discreteLevels == TRUE} and \code{colLevels} is a vector, \code{colLevels} is used as vector of breakpoints. In this case, the number of discrete colors is equal to \code{length(colLevels) - 1}.} 

\item{cutoff}{Logical. If \code{FALSE}, nodes whose \code{theme} value is beyond the range established by the vector of breakpoints are attributed the color corresponding to the lowest (or highest) value in the color scheme. If \code{TRUE}, such nodes are attributed the color \code{NaNcolor}.}

\item{colPalette}{Color palette used to display theme values. \code{colPalette} accepts both functions creating color palettes and vectors of colors. In the latter case, \code{length(colPalette)} must be greater than the number of color levels. See examples below and \code{\link{hcl.colors}}.}

\item{chooseCM}{Index of catchment to display (only effective if \code{OCN$nOutlet > 1}). It can be a logical or a numeric vector. If \code{FALSE}, all catchments are displayed. If \code{TRUE}, the catchment with largest area is displayed. If \code{chooseCM} is a subset of vector \code{1:length(OCN$nOutlet)}, only the catchment(s) identified by the indices in \code{chooseCM} are displayed.}

\item{exactDraw}{Logical. If \code{TRUE}, the real shape of OCNs is plotted. If flow crosses a boundary, the pixel that is not contiguous to its outlet is flipped.}

\item{drawNodes}{Logical. If \code{FALSE}, the theme is directly displayed on the river network. In this case, the edge departing from a given node is displayed with the color attributed to the node. If \code{TRUE}, the theme is displayed via markers at the locations of the nodes at the RN or AG level (depending on the length of \code{theme}).
In this case, \code{nanColor} can be used to define the color of the river network.}

\item{nodeType}{Only effective if \code{drawNodes == TRUE} and \code{length(theme) == OCN$RN$nNodes}. Can assume values \code{"upstream"} or \code{"downstream"}. If \code{"upstream"}, nodes are drawn at the upstream ends of the corresponding edges (i.e. at the coordinates defined by \code{OCN$AG$X}, \code{OCN$AG$Y}). If \code{"downstream"}, nodes are drawn at the downstream ends of the corresponding edges (i.e. at the coordinates defined by \code{OCN$AG$XReach}, \code{OCN$AG$YReach}).}

\item{nanColor}{Color attributed to RN or AG nodes whose theme value is \code{NA} or \code{NaN}.}

\item{riverColor}{Only effective if \code{drawNodes == TRUE}. Color used to display the OCN below the nodes.}

\item{backgroundColor}{Color used in the background of the figure. It can be either a single value, or a vector with number of components 
equal to \code{length(chooseCM)}. If \code{length(backgroundColor) == length(chooseCM)}, each color is used to identify a different catchment selected in \code{chooseCM}
(corresponding to the respective outlet). If instead \code{length(chooseCM) > 1} and \code{length(backgroundColor) == 1}, all catchments are colored with the same 
\code{backgroundColor}.}

\item{addLegend}{Logical. If \code{TRUE}, add legend to the plot. If also \code{discreteLevels = FALSE}, \code{\link{image.plot}} is used to display the legend, which appears as a colorbar; as a result, elements (e.g. node coordinates) subsequently plotted of on top of the 2D elevation map might be wrongly positioned. }

\item{min_lwd, max_lwd}{Minimum and maximum values of line width used to display the OCN (actual line width is proportional to the square root of drainage area).}

\item{...}{Further arguments to be passed to \code{plot}.}
}
\value{
No output is returned.
}
\description{
Function that draws OCNs with color of RN or AG nodes depending on an arbitrary theme.
}
\details{
This function can be used to show how a certain spatial field varies along the river network.

\emph{Default plot options.} By default, it is set \code{asp = 1}, \code{xlab = ""}, \code{ylab = ""}. If at least one between \code{xlim} and \code{ylim} is specified by the user, the default for \code{axes} is \code{TRUE}, and is \code{FALSE} if not. Specifying \code{xlim} and \code{ylim} helps zoom into a portion of the river network; however, due to the default \code{asp = 1}, the displayed region might be larger than what is expected if the ranges of \code{xlim} and \code{ylim} are different. To avoid this, set \code{asp = NA} (at the cost of producing a deformed river network).
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

\dontrun{
# 2) Display distance to outlet (at the RN level) along the main stem
# of an OCN
OCN <- aggregate_OCN(landscape_OCN(OCN_250_T)) # this takes some seconds
OCN <- paths_OCN(OCN, includePaths = TRUE) # this takes some seconds

distanceToOutlet <- OCN$RN$downstreamPathLength[,OCN$RN$outlet]
farthestNode <- which(distanceToOutlet == max(distanceToOutlet))
mainStem <- OCN$RN$downstreamPath[[farthestNode]][[OCN$RN$outlet]]
theme <- rep(NaN, OCN$RN$nNodes)
theme[mainStem] <- distanceToOutlet[mainStem]

draw_thematic_OCN(theme, OCN)
title("Distance to outlet along the main stem [pixel units]")
}

# 3) Show an OCN without a theme
OCN <- aggregate_OCN(landscape_OCN(OCN_20), thrA = 4)
draw_thematic_OCN(OCN)
draw_thematic_OCN(OCN, xlim=c(3,8), ylim=c(0,5))# zoom closer at the outlet

}
