\name{OCN_to_SSN}
\alias{OCN_to_SSN}
\title{Transform OCN into SSN object}
\usage{
OCN_to_SSN(OCN, level, obsDesign,
  predDesign = noPoints, path, importToR = FALSE)
}
\arguments{
\item{OCN}{List as produced by \code{\link{aggregate_OCN}}.}

\item{level}{Aggregation level at which the OCN is converted into a SpatialStreamNetwork object. It 
must be equal to either \code{"FD"}, \code{"RN"} or \code{"AG"}.}

\item{obsDesign}{Same as the argument of the same name in \code{\link{createSSN}}.
Note that the length of the argument of the design function must be equal to \code{OCN$N_outlet}.}

\item{predDesign}{Same as the argument of the same name in \code{\link{createSSN}}. Note that, if
a design function is specified, the length of its argument must be equal to \code{OCN$N_outlet}.}

\item{path}{Same as the argument of the same name in \code{\link{createSSN}}.}

\item{importToR}{Same as the argument of the same name in \code{\link{createSSN}}.}
}
\value{
A \code{SpatialStreamNetwork} object if \code{importToR} is \code{TRUE}, otherwise \code{NULL}.
}
\description{
Function that transforms an OCN into a SpatialStreamNetwork object. It is analogous to function 
\code{\link{createSSN}} from package \code{\link{SSN}}.
}
\details{
The generated SpatialStreamNetwork object consists of \code{OCN$N_outlet} networks. Note that an 
error is thrown if, at the selected aggregation \code{level}, at least one of these networks is degenerate
(i.e. it has less than two nodes). This is typically the case for OCNs generated with option 
\code{N_outlet = "All"}.

If \code{OCN$PeriodicBoundaries == FALSE}, nodes' locations in the SpatialStreamNetwork object are 
given by the lattice coordinates (i.e. \code{OCN$level$X}, \code{OCN$level$Y}); if 
\code{OCN$PeriodicBoundaries == TRUE}, real coordinates are used (i.e. those defined by
 \code{OCN$FD$X_draw}, \code{OCN$FD$Y_draw}, see \code{\link{landscape_OCN}}).}
\examples{
OCN <- aggregate_OCN(landscape_OCN(OCN_20), A_thr = 4)
ssn1 <- OCN_to_SSN(OCN, level = "AG", obsDesign = SSN::poissonDesign(10), 
path=paste(tempdir(),"/OCN.ssn", sep = ""), importToR = TRUE)
plot(ssn1)
}