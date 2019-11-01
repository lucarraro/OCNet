\name{create_peano}
\alias{create_peano}
\title{Create Peano network}
\usage{
create_peano(n_iter_peano, OutletPos = "NE", xllcorner = 1,
  yllcorner = 1, cellsize = 1)
}
\arguments{
\item{n_iter_peano}{Number of iteration of the Peano scheme. The resulting network will span a domain 
of size \code{2^(n_iter_peano + 1)} by \code{2^(n_iter_peano + 1)}.}

\item{OutletPos}{Corner where the outlet is located, expressed as intercardinal direction. 
Possible values are \code{"NE"}, \code{"SE"}, \code{"SW"}, \code{"NW"}.}

\item{xllcorner}{X coordinate of the lower-left pixel (expressed in planar units).}

\item{yllcorner}{Y coordinate of the lower-left pixel (expressed in planar units).}

\item{cellsize}{Size of a pixel (expressed in planar units).}
}
\value{
A list that contains the same objects as those produced by \code{\link{create_OCN}}. 
As such, it can be used as input for all other complementary functions of the package.
}
\description{
Function that creates Peano networks on a square lattice.
}
\examples{
peano <- create_peano(4)

}
