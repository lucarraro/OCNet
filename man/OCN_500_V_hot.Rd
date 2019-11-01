\docType{data}
\name{OCN_500_V_hot}
\alias{OCN_500_V_hot}
\title{Example of single-outlet OCN}
\format{A list. See \code{\link{create_OCN}} documentation for details.}
\usage{
data(OCN_500_V_hot)
}
\description{
A network built on a 500x500 lattice obtained by executing \code{set.seed(9); create_OCN(500,500,type_initialstate = "V",CoolingRate = 20,InitialNoCoolingPhase = 0.15,cellsize = 50, N_iter = 40*500*500}.
}
\keyword{datasets}
