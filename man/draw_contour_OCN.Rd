\name{draw_contour_OCN}
\alias{draw_contour_OCN}
\title{Draw Optimal Channel Network with catchment contours}
\usage{
draw_contour_OCN(OCN, A_thr_draw = 0.002 * OCN$dimX * OCN$dimY *
  OCN$cellsize^2, ExactDraw = TRUE, ColPal_river = NULL, 
  ColPal_cont = "#000000", DrawOutlets = 0, pch = 22, col = "#000000", 
  bg = "#000000")
}
\arguments{
\item{OCN}{List as produced by \code{\link{landscape_OCN}}.}

\item{A_thr_draw}{Threshold drainage area value used to display the network.}

\item{ExactDraw}{If \code{TRUE}, the real shape of OCNs is plotted. If flow crosses a boundary, 
the pixel that is not contiguous to its outlet is flipped. It is only effective if \code{OCN$PeriodicBoundaries = TRUE}}

\item{ColPal_river}{Color palette used to plot the river network(s). Default is a rearranged version of theme \code{"Dark 3"} (see \code{\link{hcl.pals}}). \code{ColPal_river} accepts both functions creating color palettes and vectors of colors (of which the first \code{OCN$N_outlet} elements are used). If a single color value is provided and \code{OCN$N_outlet > 1}, all river networks are drawn with the same color.}

\item{ColPal_cont}{Color palette used to plot the catchment contour(s). Details as in \code{ColPal_river}. Additionally, if \code{ColPal_cont = 0}, the palette specified in \code{ColPal_river} is copied. }

\item{DrawOutlets}{If equal to \code{1}, black squares are drawn at the outlets' locations behind the river; if \code{2} they are plotted on top of the river.}

\item{pch, col, bg}{Parameters controlling the way outlets are plotted (if \code{DrawOutlets = TRUE}). They correspond to the arguments of the same name of function \code{\link{points}}.}

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
# draw contour of a 20x20 single-outlet OCN
# (for single-outlet OCNs without periodic boundaries, the output
# of draw_contour_OCN is very similar to that of draw_simple_OCN)
draw_contour_OCN(landscape_OCN(OCN_20), A_thr_draw = 4)

\dontrun{
# plot real shape of multiple-outlet OCN created with periodic boundaries
OCN <- landscape_OCN(OCN_300_4out_PB) # it takes around one minute
draw_contour_OCN(OCN, DrawOutlets = 2)
}
}
