\name{draw_contour_OCN}
\alias{draw_contour_OCN}
\title{Draw Optimal Channel Network with catchment contours}
\usage{
draw_contour_OCN(OCN, A_thr_draw = 0.002 * OCN$dimX * OCN$dimY *
  OCN$cellsize^2, ExactDraw = TRUE, DrawOutlets=FALSE)
}
\arguments{
\item{OCN}{List as produced by \code{\link{landscape_OCN}}.}

\item{A_thr_draw}{Threshold drainage area value used to display the network.}

\item{ExactDraw}{If \code{TRUE}, the real shape of OCNs is plotted. If flow crosses a boundary, 
the pixel that is not contiguous to its outlet is flipped.}

\item{DrawOutlets}{If \code{TRUE}, black squares are drawn at the outlets' locations.}
}
\value{
No output is returned.
}
\description{
Function that plots real-shaped OCN and catchment contours.
}
\details{
For not too large networks (i.e. if \code{OCN$FD$Nnodes <= 40000}, corresponding to a 200x200 lattice), 
pixels whose drainage area \code{OCN$FD$A} is lower than \code{A_thr_draw} are drawn with a light grey stroke.
If \code{OCN$FD$Nnodes > 40000}, in order to speed up the execution of this function, only the network constituted 
by pixels such that \code{OCN$FD$A > A_thr_draw} is drawn.
}
\examples{
draw_contour_OCN(landscape_OCN(OCN_20))

# draw real-shaped OCN starting from an OCN with periodic boundaries:
\dontrun{
draw_contour_OCN(landscape_OCN(OCN_300_4out_PB))
}
}
