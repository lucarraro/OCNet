\docType{data}
\name{OCN_500_hot}
\alias{OCN_500_hot}
\title{Example of single-outlet OCN}
\format{A list. See \code{\link{create_OCN}} documentation for details.}
\usage{
data(OCN_500_hot)
}
\description{
A network built on a 500x500 lattice obtained by executing \code{set.seed(9); create_OCN(500, 500, outletPos = 100, outletSide = "E", coolingRate = 0.5, initialNoCoolingPhase = 0.1, cellsize = 0.5, nIter = 40*500*500)}.
}
\keyword{datasets}
