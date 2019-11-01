\name{rivergeometry_OCN}
\alias{rivergeometry_OCN}
\title{River geometry of an Optimal Channel Network}
\usage{
rivergeometry_OCN(OCN, width_max = 1, depth_max = 1,
  velocity_max = 1, exp_width = NaN, exp_depth = NaN,
  exp_velocity = NaN)
}
\arguments{
\item{OCN}{List as produced by \code{\link{aggregate_OCN}}.}

\item{width_max}{Maximum river width allowed. If \code{N_outlet = 1}, it corresponds to the width at the outlet node.}

\item{depth_max}{Maximum river depth allowed. If \code{N_outlet = 1}, it corresponds to the depth at the outlet node.}

\item{velocity_max}{Maximum water velocity allowed. If \code{N_outlet = 1}, it corresponds to the water velocity at the outlet node.}

\item{exp_width, exp_depth, exp_velocity}{Exponents for the power law relationship between river width, depth, water velocity 
and contributing area. If none of \code{exp_width}, \code{exp_depth}, \code{exp_velocity} is specified by the user, the values 
\code{exp_width = 0.5}, \code{exp_depth = 0.4}, \code{exp_depth = 0.1} proposed by Leopold and Maddock [1953] are used. 
It is possible to specify two out of these three exponents, provided that each of them lies in the range (0; 1) and their sum is lower than one. 
In this case, the missing exponent is calculated as the complement to one of the sum of the two values provided. 
If all three exponents are specified by the user, their sum must be equal to one.}
}
\value{
A list that contains all objects contained in \code{OCN}, in addition to the objects listed below. 
\item{\code{RN$width}}{Vector (of length \code{OCN$RN$Nnodes}) of river width values for every RN node.}
\item{\code{RN$depth}}{Vector (of length \code{OCN$RN$Nnodes}) of river depth values for every RN node.}
\item{\code{RN$velocity}}{Vector (of length \code{OCN$RN$Nnodes}) of water velocity values for every RN node.}
\item{\code{AG$width}}{Vector (of length \code{OCN$AG$Nnodes}) of river width values for every AG node.}
\item{\code{AG$depth}}{Vector (of length \code{OCN$AG$Nnodes}) of river depth values for every AG node.}
\item{\code{AG$velocity}}{Vector (of length \code{OCN$AG$Nnodes}) of water velocity values for every AG node.}
Finally, \code{width_max}, \code{depth_max}, \code{velocity_max}, \code{exp_width}, \code{exp_depth}, \code{exp_velocity} are added to the list.
}
\description{
Function that calculates river width, depth and water velocity by applying Leopold's scaling relationships to nodes at the RN and AG levels.
}
\details{
See also Leopold, L. B., & Maddock, T. (1953). \emph{The hydraulic geometry of stream channels and some physiographic implications}
 (Vol. 252). US Government Printing Office.
}
\examples{
OCN <- rivergeometry_OCN(aggregate_OCN(landscape_OCN(OCN_20)))
}
