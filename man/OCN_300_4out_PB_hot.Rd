\docType{data}
\name{OCN_300_4out_PB_hot}
\alias{OCN_300_4out_PB_hot}
\title{Example of multiple-outlet OCN with periodic boundaries}
\format{A list. See \code{\link{create_OCN}} documentation for details.}
\usage{
data(OCN_300_4out_PB_hot)
}
\description{
A network built on a 300x300 lattice obtained by executing \code{set.seed(5); create_OCN(300, 300, nOutlet = 4, outletSide = c("S", "N", "W", "E"), outletPos = c(1, 300, 149, 150), typeInitialState = "V", periodicBoundaries = TRUE, cellsize = 50, coolingRate = 0.5, initialNoCoolingPhase = 0.1)}.
}
\keyword{datasets}
