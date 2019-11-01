\name{landscape_OCN}
\alias{landscape_OCN}
\title{Generate 3D landscape from an Optimal Channel Network}
\usage{
landscape_OCN(OCN, Slope_0 = 1, Z_min = 0, OptimizeDZ = FALSE,
  OptimMethod = "BFGS", OptimControl = list(maxit = 100 *
  length(OCN$FD$Outlet), trace = 1), DisplayUpdates = 0)
}
\arguments{
\item{OCN}{List as produced by \code{\link{create_OCN}}.}

\item{Slope_0}{Slope of the outlet pixel (in elevation units/planar units).}

\item{Z_min}{Elevation of the lowest pixel (in elevation units).}

\item{OptimizeDZ}{If \code{TRUE}, when there are multiple catchments, minimize differences in elevation at the catchment borders 
by lifting catchments, while respecting \code{Z_min}. If \code{FALSE}, all outlet pixels have elevation equal to \code{Z_min}.}

\item{OptimMethod}{Optimization method used by function \code{\link{optim}} (only used if \code{OptimizeDZ = TRUE}).}

\item{OptimControl}{List of control parameters used by function \code{\link{optim}} (only used if \code{OptimizeDZ = TRUE}).}

\item{DisplayUpdates}{State if updates are printed on the console while \code{landscape_OCN} runs.
\describe{
\item{\code{0}}{No update is given.}
\item{\code{1}}{Concise updates are given.}
\item{\code{2}}{More extensive updates are given (this might slow down the total function runtime).}
}
Note that the display of updates during optimization of elevations (when \code{OptimizeDZ = TRUE}) is controlled by parameter \code{OptimControl$trace}.
}
}
\value{
A list that contains all objects contained in \code{OCN}, in addition to the objects listed below. 
A new sublist CM, containing variables at the catchment aggregation levels, is created.
\item{\code{FD$Slope}}{Vector (of length \code{OCN$FD$Nnodes}) of slope values (in elevation units/planar units) 
for each FD pixel, as derived by the slope/area relationship.}
\item{\code{FD$Length}}{Vector (of length \code{OCN$FD$Nnodes}) of pixel lengths. \code{OCN$FD$Length[i] = OCN$FD$cellsize} 
if flow direction in \code{i} is horizontal or vertical; \code{OCN$FD$Length[i] = OCN$FD$cellsize*sqrt(2)} if flow direction in \code{i} is diagonal.}  
\item{\code{FD$to_CM}}{Vector (of length \code{OCN$FD$Nnodes}) with catchment index values for each FD pixel. 
Example: \code{OCN$FD$to_CM[i] = j} if pixel \code{i} drains into the outlet whose location is defined by 
\code{OutletSide[j]}, \code{OutletPos[j]}.} 
\item{\code{FD$X_draw}}{When \code{PeriodicBoundaries = TRUE}, it is a vector (of length \code{OCN$FD$Nnodes}) with real X coordinate of FD pixels. 
If \code{PeriodicBoundaries = FALSE}, it is equal to \code{OCN$FD$X}.}
\item{\code{FD$Y_draw}}{When \code{PeriodicBoundaries = TRUE}, it is a vector (of length \code{OCN$FD$Nnodes}) with real Y coordinate of FD pixels. 
If \code{PeriodicBoundaries = FALSE}, it is equal to \code{OCN$FD$Y}.}
\item{\code{FD$Z}}{Vector (of length \code{OCN$FD$Nnodes}) of elevation values for each FD pixel. Values are calculated by consecutive implementation 
of the slope/area relationship along upstream paths.}
\item{\code{CM$A}}{Vector (of length \code{OCN$N_outlet}) with values of drainage area (in square planar units) for each of the catchments 
identified by the corresponding \code{OCN$FD$Outlet}.}
\item{\code{CM$W}}{Adjacency matrix (\code{OCN$N_outlet} by \code{OCN$N_outlet}) at the catchment level. Two catchments are connected 
if they share a border. Note that this is not a flow connection. Unlike the adjacency matrices at levels FD, RN, AG, this matrix is symmetric. 
It is a \code{\link{spam}} object.}
\item{\code{CM$X_contour} (\code{CM$Y_contour})}{List with number of objects equal to \code{OCN$FD$N_outlet}. Each object \code{i} is a list 
with X (Y) coordinates of the contour of catchment \code{i} for use in plots with \code{ExactDraw = FALSE} (see functions \code{\link{draw_contour_OCN}}, 
\code{\link{draw_thematic_OCN}}). If catchment \code{i} is constituted by regions that are only connected through a diagonal flow direction, 
\code{CM$X_contour[[i]]} (\code{CM$Y_contour[[i]]}) contains as many objects as the number of regions into which catchment \code{i} is split.} 
\item{\code{CM$X_contour_draw} (\code{CM$Y_contour_draw})}{List with number of objects equal to \code{OCN$FD$N_outlet}. Each object \code{i} is a list 
with X (Y) coordinates of the contour of catchment \code{i} for use in plots with \code{ExactDraw = TRUE} (see functions \code{\link{draw_contour_OCN}}, 
\code{\link{draw_thematic_OCN}}). If catchment \code{i} is constituted by regions that are only connected through a diagonal flow direction, 
\code{CM$X_contour_draw[[i]]} (\code{CM$Y_contour_draw[[i]]}) contains as many objects as the number of regions into which catchment \code{i} is split.}
\item{\code{OptList}}{List of output parameters produced by the optimization function \code{\link{optim}} (only present if \code{OptimizeDZ = TRUE}).}
Finally, \code{Slope_0} and \code{Z_min} are passed to the list as they were included in the input.
}
\description{
Function that calculates the elevation field generated by the OCN and the partition of the domain into different catchments.
}
\details{
The function features an algorithm (which can be activated via the optional input \code{OptimizeDZ}) that, given the network 
configuration and a \code{Slope_0} value, finds the elevation of \code{OCN$N_outlet - 1} outlets relative to the elevation of the first 
outlet in vectors \code{OutletSide}, \code{OutletPos} such that the sum of the absolute differences elevation of neighboring pixels
belonging to different catchments is minimized. Such objective function is minimized by means of function \code{\link{optim}}. 
The absolute elevation of the outlet pixels (and, consequently, of the whole lattice) is finally attributed by imposing 
\code{OCN$FD$Z >= Z_min}. Note that, due to the high dimensionality of the problem, convergence of the 
optimization algorithm is not guaranteed for large \code{OCN$N_outlet} (say, \code{OCN$N_outlet > 10}).
}
\examples{
OCN2 <- landscape_OCN(OCN_20)
}
