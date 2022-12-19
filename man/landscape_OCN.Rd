\name{landscape_OCN}
\alias{landscape_OCN}
\title{Generate 3D landscape from an Optimal Channel Network}
\usage{
landscape_OCN(OCN, slope0 = 1, zMin = 0, optimizeDZ = FALSE,
  optimMethod = "BFGS", optimControl = list(maxit = 100 *
  length(OCN$FD$outlet), trace = 1), displayUpdates = 0)
}
\arguments{
\item{OCN}{List as produced by \code{\link{create_OCN}}.}

\item{slope0}{slope of the outlet pixel (in elevation units/planar units).}

\item{zMin}{Elevation of the lowest pixel (in elevation units).}

\item{optimizeDZ}{If \code{TRUE}, when there are multiple catchments, minimize differences in elevation at the catchment borders 
by lifting catchments, while respecting \code{zMin}. If \code{FALSE}, all outlet pixels have elevation equal to \code{zMin}. 
This option is not effective for OCNs generated via \code{\link{create_general_contour_OCN}}.}

\item{optimMethod}{Optimization method used by function \code{\link{optim}} (only used if \code{optimizeDZ = TRUE}).}

\item{optimControl}{List of control parameters used by function \code{\link{optim}} (only used if \code{optimizeDZ = TRUE}).}

\item{displayUpdates}{State if updates are printed on the console while \code{landscape_OCN} runs.
\describe{
\item{\code{0}}{No update is given.}
\item{\code{1}}{Concise updates are given.}
\item{\code{2}}{More extensive updates are given (this might slow down the total function runtime).}
}
Note that the display of updates during optimization of elevations (when \code{optimizeDZ = TRUE}) is controlled by parameter \code{optimControl$trace}.
}
}
\value{
A list that contains all objects contained in \code{OCN}, in addition to the objects listed below. 
A new sublist CM, containing variables at the catchment aggregation levels, is created.
\item{\code{FD$slope}}{Vector (of length \code{OCN$FD$nNodes}) of slope values (in elevation units/planar units) 
for each FD pixel, as derived by the slope/area relationship.}
\item{\code{FD$leng}}{Vector (of length \code{OCN$FD$nNodes}) of pixel lengths. \code{OCN$FD$leng[i] = OCN$FD$cellsize} 
if flow direction in \code{i} is horizontal or vertical; \code{OCN$FD$leng[i] = OCN$FD$cellsize*sqrt(2)} if flow direction in \code{i} is diagonal.}  
\item{\code{FD$toCM}}{Vector (of length \code{OCN$FD$nNodes}) with catchment index values for each FD pixel. 
Example: \code{OCN$FD$toCM[i] = j} if pixel \code{i} drains into the outlet whose location is defined by 
\code{outletSide[j]}, \code{outletPos[j]}.} 
\item{\code{FD$XDraw}}{When \code{periodicBoundaries = TRUE}, it is a vector (of length \code{OCN$FD$nNodes}) with real X coordinate of FD pixels. 
If \code{periodicBoundaries = FALSE}, it is equal to \code{OCN$FD$X}.}
\item{\code{FD$YDraw}}{When \code{periodicBoundaries = TRUE}, it is a vector (of length \code{OCN$FD$nNodes}) with real Y coordinate of FD pixels. 
If \code{periodicBoundaries = FALSE}, it is equal to \code{OCN$FD$Y}.}
\item{\code{FD$Z}}{Vector (of length \code{OCN$FD$nNodes}) of elevation values for each FD pixel. Values are calculated by consecutive implementation 
of the slope/area relationship along upstream paths.}
\item{\code{CM$A}}{Vector (of length \code{OCN$nOutlet}) with values of drainage area (in square planar units) for each of the catchments 
identified by the corresponding \code{OCN$FD$outlet}.}
\item{\code{CM$W}}{Adjacency matrix (\code{OCN$nOutlet} by \code{OCN$nOutlet}) at the catchment level. Two catchments are connected 
if they share a border. Note that this is not a flow connection. Unlike the adjacency matrices at levels FD, RN, AG, this matrix is symmetric. 
It is a \code{\link{spam}} object.}
\item{\code{CM$XContour} (\code{CM$Y_contour})}{List with number of objects equal to \code{OCN$FD$nOutlet}. Each object \code{i} is a list 
with X (Y) coordinates of the contour of catchment \code{i} for use in plots with \code{exactDraw = FALSE} (see functions \code{\link{draw_contour_OCN}}, 
\code{\link{draw_thematic_OCN}}). If catchment \code{i} is constituted by regions that are only connected through a diagonal flow direction, 
\code{CM$XContour[[i]]} (\code{CM$Y_contour[[i]]}) contains as many objects as the number of regions into which catchment \code{i} is split.} 
\item{\code{CM$XContourDraw} (\code{CM$YContourDraw})}{List with number of objects equal to \code{OCN$FD$nOutlet}. Each object \code{i} is a list 
with X (Y) coordinates of the contour of catchment \code{i} for use in plots with \code{exactDraw = TRUE} (see functions \code{\link{draw_contour_OCN}}, 
\code{\link{draw_thematic_OCN}}). If catchment \code{i} is constituted by regions that are only connected through a diagonal flow direction, 
\code{CM$XContourDraw[[i]]} (\code{CM$YContourDraw[[i]]}) contains as many objects as the number of regions into which catchment \code{i} is split.}
\item{\code{OptList}}{List of output parameters produced by the optimization function \code{\link{optim}} (only present if \code{optimizeDZ = TRUE}).}
Finally, \code{slope0} and \code{zMin} are passed to the list as they were included in the input.
}
\description{
Function that calculates the elevation field generated by the OCN and the partition of the domain into different catchments.
}
\details{
The function features an algorithm (which can be activated via the optional input \code{optimizeDZ}) that, given the network 
configuration and a \code{slope0} value, finds the elevation of \code{OCN$nOutlet - 1} outlets relative to the elevation of the first 
outlet in vectors \code{outletSide}, \code{outletPos} such that the sum of the absolute differences elevation of neighboring pixels
belonging to different catchments is minimized. Such objective function is minimized by means of function \code{\link{optim}}. 
The absolute elevation of the outlet pixels (and, consequently, of the whole lattice) is finally attributed by imposing 
\code{OCN$FD$Z >= zMin}. Note that, due to the high dimensionality of the problem, convergence of the 
optimization algorithm is not guaranteed for large \code{OCN$nOutlet} (say, \code{OCN$nOutlet > 10}).
}
\examples{
# 1) draw 2D elevation map of a 20x20 OCN with default options
OCN2 <- landscape_OCN(OCN_20)

\dontrun{
# 2) generate a 100x50 OCN; assume that the pixel resolution is 200 m
# (the total catchment area is 20 km2)
set.seed(1)
OCN <- create_OCN(100, 50, cellsize = 200, 
		displayUpdates = 0) # this takes about 40 s
# use landscape_OCN to derive the 3D landscape subsumed by the OCN
# by assuming that the elevation and slope at the outlet are 200 m 
# and 0.0075, respectively
OCN <- landscape_OCN(OCN, zMin = 200, slope0 = 0.0075)
# draw 2D and 3D representations of the landscape
draw_elev2D_OCN(OCN)
draw_elev3D_OCN(OCN)
draw_elev3Drgl_OCN(OCN)
}

\dontrun{
# 3) generate a 100x50 OCN with 4 outlets
set.seed(1)
OCN <- create_OCN(100, 50, cellsize = 200, 
		nOutlet = 4, displayUpdates = 0) # this takes about 40 s
# use landscape_OCN and optimize elevation of outlets	
OCN <- landscape_OCN(OCN, slope0 = 0.0075, 
		optimizeDZ = TRUE)
# display elevation of outlets and 2D elevation map
OCN$FD$Z[OCN$FD$outlet]	
draw_elev2D_OCN(OCN)	
}
}
