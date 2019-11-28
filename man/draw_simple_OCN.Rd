\name{draw_simple_OCN}
\alias{draw_simple_OCN}
\title{Draw an Optimal Channel Network}
\usage{
draw_simple_OCN(OCN, A_thr_draw = 0.002 * OCN$dimX * OCN$dimY *
  OCN$cellsize^2, plot_title = "", river_color = "#0066FF", EasyDraw = NULL)
}
\arguments{
\item{OCN}{List as produced by \code{\link{create_OCN}}.}

\item{A_thr_draw}{Threshold drainage area value used to display the network.}

\item{plot_title}{String that will appear as title of the plot.}

\item{river_color}{Color used to plot the river.}

\item{EasyDraw}{Logical variable. If \code{TRUE}, the whole network is displayed, and pixels with drainage area lower than \code{A_thr_draw} are displayed in light gray. If \code{FALSE}, only pixels with drainage area greater or equal to \code{A_thr_draw} are displayed. Default is \code{FALSE} if \code{OCN$Nnodes <= 40000}, and \code{TRUE} otherwise. Note that setting \code{EasyDraw = FALSE} for large networks might slow down the process considerably.}
}
\value{
No output is returned.
}
\description{
Function that plots the non-aggregated OCN as calculated by \code{create_OCN}.
}
\examples{
# draw OCN with default settings
draw_simple_OCN(OCN_250)
# same as above, but with decreased A_thr_draw
draw_simple_OCN(OCN_250, A_thr_draw = 0.001 * OCN_250$dimX * OCN_250$dimY)
# same as the first example, but include the portion of network  
# with drainage area lower than A_thr_draw 
\dontrun{
draw_simple_OCN(OCN_250, EasyDraw = FALSE) # this will take some seconds
}
}
