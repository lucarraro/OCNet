\name{OCN_to_igraph}
\alias{OCN_to_igraph}
\title{Transform OCN into igraph object}
\usage{
OCN_to_igraph(OCN, level)
}
\arguments{
\item{OCN}{A \code{river} object as produced by \code{\link{aggregate_OCN}}.}

\item{level}{Aggregation level at which the OCN is converted into an igraph object. It 
must be equal to either \code{"FD"}, \code{"RN"} or \code{"AG"}.}
}
\value{
An \code{igraph} object.
}
\description{
Function that transforms an OCN into an igraph object. 
}
\examples{
# 1) transform a 20x20 OCN, at the AG level, into a graph object 
OCN <- aggregate_OCN(landscape_OCN(OCN_20), thrA = 4)
g <- OCN_to_igraph(OCN, level = "AG")
plot(g, layout = matrix(c(OCN$AG$X,OCN$AG$Y), ncol = 2, nrow = OCN$AG$nNodes))
}
