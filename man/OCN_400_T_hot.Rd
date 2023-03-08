\docType{data}
\name{OCN_400_T_hot}
\alias{OCN_400_T_hot}
\title{Example of single-outlet OCN}
\format{A \code{river} object. See \code{\link{create_OCN}} documentation for details.}
\usage{
data(OCN_400_T_hot)
}
\description{
A network built on a 400x400 lattice obtained by executing \code{set.seed(7); create_OCN(400, 400, typeInitialState = "T", coolingRate = 0.5, initialNoCoolingPhase = 0.1, cellsize = 50)}.
}
\keyword{datasets}
