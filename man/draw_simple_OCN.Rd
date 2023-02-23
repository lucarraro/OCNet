\name{draw_simple_OCN}
\alias{draw_simple_OCN}
\title{Draw an Optimal Channel Network}
\usage{
draw_simple_OCN(OCN, thrADraw = 0.002 * OCN$FD$nNodes *
  OCN$cellsize^2, riverColor = "#0066FF", easyDraw = NULL, 
  min_lwd = 0.5,  max_lwd = 5)
}
\arguments{
\item{OCN}{List as produced by \code{\link{create_OCN}}.}

\item{thrADraw}{Threshold drainage area value used to display the network.}

\item{riverColor}{Color used to plot the river.}

\item{easyDraw}{Logical. If \code{TRUE}, the whole network is displayed, and pixels with drainage area lower than \code{thrADraw} are displayed in light gray. If \code{FALSE}, only pixels with drainage area greater or equal to \code{thrADraw} are displayed. Default is \code{FALSE} if \code{OCN$nNodes <= 40000}, and \code{TRUE} otherwise. Note that setting \code{easyDraw = FALSE} for large networks might slow down the process considerably.}

\item{min_lwd, max_lwd}{Minimum and maximum values of line width used to display the OCN (actual line width is proportional to the square root of drainage area).}
}
\value{
No output is returned.
}
\description{
Function that plots the non-aggregated OCN as calculated by \code{create_OCN}.
}
\examples{
# 1a) draw OCN with default settings
draw_simple_OCN(OCN_250_T)
# 1b) same as above, but with decreased thrADraw
draw_simple_OCN(OCN_250_T, thrADraw = 0.001 * OCN_250_T$dimX * OCN_250_T$dimY)
\donttest{
# 1c) same as the first example, but include the portion of network  
# with drainage area lower than thrADraw 
draw_simple_OCN(OCN_250_T, easyDraw = FALSE) # this will take some seconds
}
}
