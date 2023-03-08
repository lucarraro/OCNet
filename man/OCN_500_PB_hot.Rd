\docType{data}
\name{OCN_500_PB_hot}
\alias{OCN_500_PB_hot}
\title{Example of single-outlet OCN with periodic boundaries}
\format{A \code{river} object. See \code{\link{create_OCN}} documentation for details.}
\usage{
data(OCN_500_PB_hot)
}
\description{
A network built on a 500x500 lattice obtained by executing \code{set.seed(10); create_OCN(500, 500, periodicBoundaries = TRUE, coolingRate = 0.5, initialNoCoolingPhase = 0.1, cellsize = 50)}.
}
\keyword{datasets}
