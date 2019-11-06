\name{find_area_threshold_OCN}
\alias{find_area_threshold_OCN}
\title{Find relationship between number of nodes and threshold area in an OCN}
\usage{
find_area_threshold_OCN(OCN, thr_values = seq(OCN$cellsize,
  OCN$cellsize * max(OCN$FD$A), OCN$cellsize), MaxReachLength = Inf,
  StreamOrderType = "Strahler", DisplayUpdates = 0)
}
\arguments{
\item{OCN}{List as produced by \code{\link{landscape_OCN}}}

\item{thr_values}{Vector of values of threshold pixel area (in squared planar units) for which the respective number of nodes is computed.}

\item{MaxReachLength}{Maximum reach length allowed (in planar units). If the path length between a channel head and the downstream confluence 
is higher than \code{MaxReachLength}, the reach starting from the channel head will have a length up to \code{MaxReachLength}, while the next downstream 
pixel is considered as a new channel head, from which a new reach departs.}

\item{StreamOrderType}{If \code{"Strahler"}, Strahler stream order is computed; if \code{"Shreve"}, Shreve stream order is computed.}

\item{DisplayUpdates}{If \code{1}, progress updates are printed in the console while the function is running. If \code{0}, no updates are printed.}
}
\value{
A list whose objcets are listed below.
\item{\code{thr_values}}{Copy of the homonymous input vector.}
\item{\code{Nnodes}}{Vector (of the same length as \code{thr_values}) of number of nodes resulting from the aggregation process 
with a threshold area values specified by \code{thr_values}.}
\item{\code{DrainageDensity}}{Vector (of the same length as \code{thr_values}) of values of drainage density of the river network 
resulting from the aggregation process with a threshold area values specified by \code{thr_values}. Drainage density is calculated 
as total length of the river network divided by area of the lattice. It is expressed in planar units^(-1).}
\item{\code{StreamOrder}}{Vector (of the same length as \code{thr_values}) of values of maximum stream order attained by the river network, 
resulting from the aggregation process with a threshold area values specified by \code{thr_values}.}
}
\description{
Function that calculates relationship between threshold area and number of nodes at AG level for a given OCN. 
It can be used prior to application of \code{\link{aggregate_OCN}} in order to derive the drainage area threshold 
that corresponds to the desired number of nodes of the aggregated network.
}
\examples{
thr <- find_area_threshold_OCN(landscape_OCN(OCN_20))
}
