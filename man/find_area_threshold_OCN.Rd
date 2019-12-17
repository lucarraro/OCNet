\name{find_area_threshold_OCN}
\alias{find_area_threshold_OCN}
\title{Find relationship between number of nodes and threshold area in an OCN}
\usage{
find_area_threshold_OCN(OCN, thrValues = seq(OCN$cellsize^2,
  min(OCN$CM$A), OCN$cellsize^2), maxReachLength = Inf,
  streamOrderType = "Strahler", displayUpdates = 0)
}
\arguments{
\item{OCN}{List as produced by \code{\link{landscape_OCN}}}

\item{thrValues}{Vector of values of threshold drainage area (in squared planar units) for which the respective number of nodes at the RN and AG levels are computed. Note that it must be \code{max(thrValues) <= min(OCN$CM$A)}, otherwise the catchment(s) with area lower than \code{max(thrValues)} degenerate to a network with zero nodes at the RN/AG level.}

\item{maxReachLength}{Maximum reach length allowed (in planar units). If the path length between a channel head and the downstream confluence 
is higher than \code{maxReachLength}, the reach starting from the channel head will have a length up to \code{maxReachLength}, while the next downstream 
pixel is considered as a new channel head, from which a new reach departs.}

\item{streamOrderType}{If \code{"Strahler"}, Strahler stream order is computed; if \code{"Shreve"}, Shreve stream order is computed.}

\item{displayUpdates}{If \code{1}, progress updates are printed in the console while the function is running. If \code{0}, no updates are printed.}
}
\value{
A list whose objects are listed below.
\item{\code{thrValues}}{Copy of the input vector with the same name.}
\item{\code{nNodesRN}}{Vector (of the same length as \code{thrValues}) of number of nodes at the RN level resulting from the aggregation process 
with a threshold area values specified by \code{thrValues}.}
\item{\code{nNodesAG}}{Vector (of the same length as \code{thrValues}) of number of nodes at the AG level resulting from the aggregation process 
with a threshold area values specified by \code{thrValues}.}
\item{\code{drainageDensity}}{Vector (of the same length as \code{thrValues}) of values of drainage density of the river network 
resulting from the aggregation process with a threshold area values specified by \code{thrValues}. Drainage density is calculated 
as total length of the river network divided by area of the lattice. It is expressed in planar units^(-1).}
\item{\code{streamOrder}}{Vector (of the same length as \code{thrValues}) of values of maximum stream order attained by the river network, 
resulting from the aggregation process with a threshold area values specified by \code{thrValues}.}
}
\description{
Function that calculates relationship between threshold area and number of nodes at RN and AG level for a given OCN. 
It can be used prior to application of \code{\link{aggregate_OCN}} in order to derive the drainage area threshold 
that corresponds to the desired number of nodes of the aggregated network.
It is intended for use with single outlet OCNs, although its use with multiple outlet OCNs is allowed (provided that \code{max(thrValues) <= min(OCN$CM$A)}).
}
\examples{
# 1) derive relationship between threshold area and number of nodes
OCN <- landscape_OCN(OCN_20)
thr <- find_area_threshold_OCN(OCN)
# log-log plot of number of nodes at the AG level versus 
# relative threshold area (as fraction of total drainage area) 
old.par <- par(no.readonly = TRUE)
par(mai = c(1,1,1,1))  
plot(thr$thrValues[thr$nNodesAG > 0]/OCN$CM$A,
		thr$nNodesAG[thr$nNodesAG > 0], log = "xy",
		xlab = "Relative area threshold", ylab = "Number of AG nodes")
par(old.par)
}
