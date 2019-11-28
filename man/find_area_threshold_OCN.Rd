\name{find_area_threshold_OCN}
\alias{find_area_threshold_OCN}
\title{Find relationship between number of nodes and threshold area in an OCN}
\usage{
find_area_threshold_OCN(OCN, thr_values = seq(OCN$cellsize^2,
  min(OCN$CM$A), OCN$cellsize^2), MaxReachLength = Inf,
  StreamOrderType = "Strahler", DisplayUpdates = 0)
}
\arguments{
\item{OCN}{List as produced by \code{\link{landscape_OCN}}}

\item{thr_values}{Vector of values of threshold drainage area (in squared planar units) for which the respective number of nodes at the RN and AG levels are computed. Note that it must be \code{max(thr_values) <= min(OCN$CM$A)}, otherwise the catchment(s) with area lower than \code{max(thr_values)} degenerate to a network with zero nodes at the RN/AG level.}

\item{MaxReachLength}{Maximum reach length allowed (in planar units). If the path length between a channel head and the downstream confluence 
is higher than \code{MaxReachLength}, the reach starting from the channel head will have a length up to \code{MaxReachLength}, while the next downstream 
pixel is considered as a new channel head, from which a new reach departs.}

\item{StreamOrderType}{If \code{"Strahler"}, Strahler stream order is computed; if \code{"Shreve"}, Shreve stream order is computed.}

\item{DisplayUpdates}{If \code{1}, progress updates are printed in the console while the function is running. If \code{0}, no updates are printed.}
}
\value{
A list whose objects are listed below.
\item{\code{thr_values}}{Copy of the homonymous input vector.}
\item{\code{Nnodes_RN}}{Vector (of the same length as \code{thr_values}) of number of nodes at the RN level resulting from the aggregation process 
with a threshold area values specified by \code{thr_values}.}
\item{\code{Nnodes_AG}}{Vector (of the same length as \code{thr_values}) of number of nodes at the AG level resulting from the aggregation process 
with a threshold area values specified by \code{thr_values}.}
\item{\code{DrainageDensity}}{Vector (of the same length as \code{thr_values}) of values of drainage density of the river network 
resulting from the aggregation process with a threshold area values specified by \code{thr_values}. Drainage density is calculated 
as total length of the river network divided by area of the lattice. It is expressed in planar units^(-1).}
\item{\code{StreamOrder}}{Vector (of the same length as \code{thr_values}) of values of maximum stream order attained by the river network, 
resulting from the aggregation process with a threshold area values specified by \code{thr_values}.}
}
\description{
Function that calculates relationship between threshold area and number of nodes at RN and AG level for a given OCN. 
It can be used prior to application of \code{\link{aggregate_OCN}} in order to derive the drainage area threshold 
that corresponds to the desired number of nodes of the aggregated network.
It is intended for use with single outlet OCNs, although its use with multiple outlet OCNs is allowed (provided that \code{max(thr_values) <= min(OCN$CM$A)}).
}
\examples{
# derive relationship between threshold area and number of nodes
OCN <- landscape_OCN(OCN_20)
thr <- find_area_threshold_OCN(OCN)
# log-log plot of number of nodes at the AG level versus 
# relative threshold area (as fraction of total drainage area) 
par(mai=c(1,1,1,1))  
plot(thr$thr_values[thr$Nnodes_AG > 0]/OCN$CM$A,
		thr$Nnodes_AG[thr$Nnodes_AG > 0], log = "xy",
		xlab = "Relative area threshold", ylab = "Number of AG nodes")
}
