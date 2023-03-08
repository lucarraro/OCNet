\name{rivergeometry_OCN}
\alias{rivergeometry_OCN}
\title{River geometry of an Optimal Channel Network}
\usage{
rivergeometry_OCN(OCN, widthMax = 1, depthMax = 1,
  velocityMax = 1, expWidth = NaN, expDepth = NaN,
  expVelocity = NaN)
}
\arguments{
\item{OCN}{A \code{river} object as produced by \code{\link{aggregate_OCN}}.}

\item{widthMax}{Maximum river width allowed. If \code{nOutlet = 1}, it corresponds to the width at the outlet node.}

\item{depthMax}{Maximum river depth allowed. If \code{nOutlet = 1}, it corresponds to the depth at the outlet node.}

\item{velocityMax}{Maximum water velocity allowed. If \code{nOutlet = 1}, it corresponds to the water velocity at the outlet node.}

\item{expWidth, expDepth, expVelocity}{Exponents for the power law relationship between river width, depth, water velocity 
and contributing area. If none of \code{expWidth}, \code{expDepth}, \code{expVelocity} is specified by the user, the values 
\code{expWidth = 0.5}, \code{expDepth = 0.4}, \code{expDepth = 0.1} proposed by Leopold and Maddock [1953] are used. 
It is possible to specify two out of these three exponents, provided that each of them lies in the range (0; 1) and their sum is lower than one. 
In this case, the missing exponent is calculated as the complement to one of the sum of the two values provided. 
If all three exponents are specified by the user, their sum must be equal to one.}
}
\value{
AA \code{river} object that contains all objects contained in \code{OCN}, in addition to the objects listed below. 
\item{\code{RN$width}}{Vector (of length \code{OCN$RN$nNodes}) of river width values for every RN node.}
\item{\code{RN$depth}}{Vector (of length \code{OCN$RN$nNodes}) of river depth values for every RN node.}
\item{\code{RN$velocity}}{Vector (of length \code{OCN$RN$nNodes}) of water velocity values for every RN node.}
\item{\code{AG$width}}{Vector (of length \code{OCN$AG$nNodes}) of river width values for every AG node.}
\item{\code{AG$depth}}{Vector (of length \code{OCN$AG$nNodes}) of river depth values for every AG node.}
\item{\code{AG$velocity}}{Vector (of length \code{OCN$AG$nNodes}) of water velocity values for every AG node.}
Finally, \code{widthMax}, \code{depthMax}, \code{velocityMax}, \code{expWidth}, \code{expDepth}, \code{expVelocity} are added to the list.
}
\description{
Function that calculates river width, depth and water velocity by applying Leopold's scaling relationships to nodes at the RN and AG levels.
}
\details{
The values of contributing area used to evaluate river geometry at the AG level are equal to \code{0.5*(OCN$AG$A + OCN$AG$AReach)}. See also \code{\link{aggregate_OCN}}.

See also Leopold, L. B., & Maddock, T. (1953). \emph{The hydraulic geometry of stream channels and some physiographic implications}
 (Vol. 252). US Government Printing Office.
}
\examples{
# 1) Compute river geometry of a 20x20 OCN with default options
# and display river width at the RN level
OCN <- rivergeometry_OCN(aggregate_OCN(landscape_OCN(OCN_20)))
draw_thematic_OCN(OCN$RN$width,OCN)
}
