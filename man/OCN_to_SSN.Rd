\name{OCN_to_SSN}
\alias{OCN_to_SSN}
\title{Transform OCN into SSN object}
\usage{
OCN_to_SSN(OCN, level, obsDesign, obsSites, predDesign, 
  predSites, path, randomAllocation = FALSE, importToR = FALSE)
}
\arguments{
\item{OCN}{A \code{river} object as produced by \code{\link{aggregate_OCN}}.}

\item{level}{Aggregation level at which the OCN is converted into a SpatialStreamNetwork object. It 
must be equal to either \code{"FD"}, \code{"RN"} or \code{"AG"}.}

\item{obsDesign}{Same as the argument of the same name in \code{\link{createSSN}}.
Note that the length of the argument of the design function must be equal to \code{OCN$nOutlet}.
Only one between \code{obsDesign} and \code{obsSites} must be specified.} 

\item{obsSites}{Identifier of nodes at the aggregation \code{level} where observation sites are 
to be located. Only one between \code{obsDesign} and \code{obsSites} must be specified.}

\item{predDesign}{Same as the argument of the same name in \code{\link{createSSN}}. Note that, if
a design function is specified, the length of its argument must be equal to \code{OCN$nOutlet}.}

\item{predSites}{Identifier of nodes at the aggregation \code{level} where prediction sites are 
to be located. Both \code{predDesign} and \code{predSites} are optional, but they cannot be both specified.}

\item{path}{Same as the argument of the same name in \code{\link{createSSN}}.}

\item{randomAllocation}{Logical. If \code{FALSE}, observation and prediction nodes (when \code{obsSites} and
\code{predSites} are specified, respectively) are located at the exact position of the corresponding nodes.
If \code{TRUE}, they are located at a random position along the edge departing from the selected nodes.}

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
The generated SpatialStreamNetwork object consists of \code{OCN$nOutlet} networks. Note that an 
error is thrown if, at the selected aggregation \code{level}, at least one of these networks is degenerate
(i.e. it has less than two nodes). This is typically the case for OCNs generated with option 
\code{nOutlet = "All"}.

If \code{OCN$periodicBoundaries == FALSE}, nodes' locations in the SpatialStreamNetwork object are 
given by the lattice coordinates (i.e. \code{OCN$level$X}, \code{OCN$level$Y}); if 
\code{OCN$periodicBoundaries == TRUE}, real coordinates are used (i.e. those defined by
 \code{OCN$FD$Xdraw}, \code{OCN$FD$Ydraw}, see \code{\link{landscape_OCN}}).}
\examples{
# transform a 20x20 single-outlet OCN (aggregated at the AG level)
# into a SSN object and plot it
OCN <- aggregate_OCN(landscape_OCN(OCN_20), thrA = 4)
ssn1 <- OCN_to_SSN(OCN, "AG", obsDesign = SSN::poissonDesign(10), 
path=paste(tempdir(),"/OCN.ssn", sep = ""), importToR = TRUE)
plot(ssn1)

\donttest{
# 1) create a 50x50 OCN with two outlets and periodic boundaries; 
set.seed(1)
OCN <- create_OCN(50, 50, nOutlet = 2, outletSide = c("S", "N"),
		outletPos = c(1, 50), periodicBoundaries = TRUE)
# aggregate the OCN;
OCN <- aggregate_OCN(landscape_OCN(OCN))
# transform it into a SSN object aggregated at the RN level;
ssn2 <- OCN_to_SSN(OCN, "RN", obsDesign = SSN::binomialDesign(c(10, 10)),
		path = paste(tempdir(),"/OCN2.ssn", sep = ""), importToR = TRUE)
# and plot the SSN object; it is constituted by two networks,
# and nodes' coordinates are the "real" ones  
old.par <- par(no.readonly = TRUE)
par(mai = c(1, 1, 1, 1))
plot(ssn2, xlab = "X", ylab = "Y")
par(old.par)


# 2) An OCN with 3 outlets transformed into a SSN, with obsSites
# attributed at given locations
set.seed(2)
OCN <- aggregate_OCN(landscape_OCN(create_OCN(20, 20, nOutlet = 3)), thrA = 4)
obsSites <- c(1, 7, 12, 15, 21, 28, 34, 41, 50)
ssn3 <- OCN_to_SSN(OCN, "AG", obsSites=obsSites, 
                    path = paste(tempdir(),"/OCN3.ssn", sep = ""), importToR = TRUE)
draw_thematic_OCN(OCN)
points(OCN$AG$X[obsSites], OCN$AG$Y[obsSites])
plot(ssn3)
}
}
