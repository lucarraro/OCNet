\name{create_peano}
\alias{create_peano}
\title{Create Peano network}
\usage{
create_peano(nIterPeano, outletPos = "NE", xllcorner = 1,
  yllcorner = 1, cellsize = 1)
}
\arguments{
\item{nIterPeano}{Number of iteration of the Peano scheme. The resulting network will span a domain 
of size \code{2^(nIterPeano + 1)} by \code{2^(nIterPeano + 1)}.}

\item{outletPos}{Corner where the outlet is located, expressed as intercardinal direction. 
Possible values are \code{"NE"}, \code{"SE"}, \code{"SW"}, \code{"NW"}.}

\item{xllcorner}{X coordinate of the lower-left pixel (expressed in planar units).}

\item{yllcorner}{Y coordinate of the lower-left pixel (expressed in planar units).}

\item{cellsize}{Size of a pixel (expressed in planar units).}
}
\value{
A \code{river} object that contains the same objects as those produced by \code{\link{create_OCN}}. 
As such, it can be used as input for all other complementary functions of the package.
}
\description{
Function that creates Peano networks on a square lattice.
}
\examples{
# 1) create a peano network in a 32x32 square,
# use landscape_OCN, aggregate_OCN functions,
# and display subcatchment map and map of drainage area
peano <- create_peano(4)
peano <- aggregate_OCN(landscape_OCN(peano), thrA = 4)
old.par <- par(no.readonly = TRUE)
par(mfrow=c(1,3))
draw_simple_OCN(peano)
title("Peano network")
draw_subcatchments_OCN(peano)
title("Subcatchments")
draw_thematic_OCN(peano$RN$A, peano)
title("Drainage area at RN level")
par(old.par)
}
