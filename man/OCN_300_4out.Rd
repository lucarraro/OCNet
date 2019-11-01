\docType{data}
\name{OCN_300_4out}
\alias{OCN_300_4out}
\title{Example of multiple-outlet OCN}
\format{A list. See \code{\link{create_OCN}} documentation for details.}
\usage{
data(OCN_300_4out)
}
\description{
A network built on a 300x300 lattice obtained by executing \code{set.seed(5); create_OCN(300,300,N_outlet = 4,OutletSide = c("S","N","W","E"),OutletPos = c(1,300,149,150),type_initialstate = "V",cellsize = 50)}.
}
\keyword{datasets}
