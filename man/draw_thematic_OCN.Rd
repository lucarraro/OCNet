\name{draw_thematic_OCN}
\alias{draw_thematic_OCN}
\title{Draw thematic map on an Optimal Channel Network}
\usage{
draw_thematic_OCN(theme, OCN,
  minval = min(theme[!(is.nan(theme))]),
  maxval = max(theme[!(is.nan(theme))]), plot_title = "",
  ColPalette = colorRampPalette(c("yellow","red","black")),
  ColLevels = 1000, DiscreteLevels = FALSE, ExactDraw = FALSE,
  DrawNodes = FALSE, Node_cex = 2, Node_pch = 21,
  NaNColor = "#0099FF", BackgroundColor = NULL,
  add_colorbar=TRUE)
}
\arguments{
\item{theme}{Vector (of length \code{OCN$AG$Nnodes} or \code{OCN$RN$Nnodes}) expressing the spatial field of interest. 
The vector can contain \code{NA} and \code{NaN} values to identify RN or AG nodes where the theme is not defined.}

\item{OCN}{List as produced by \code{\link{aggregate_OCN}}.}

\item{minval}{Lower bound of the color scale.}

\item{maxval}{Upper bound of the color scale.}

\item{plot_title}{String that will appear as title of the plot.}

\item{ColPalette}{Color palette used to display theme values. \code{ColPalette} accepts both functions creating color palettes and vectors of colors (of which the first \code{ColLevels} elements are used). See \code{\link{Palettes}} and examples below.}

\item{ColLevels}{Number of colors in the palette.}

\item{DiscreteLevels}{If \code{TRUE}, \code{ColLevels} is set equal to the number of unique values in theme. 
In this way, the colorbar displays only the colors that are actually used in the map.}

\item{ExactDraw}{If \code{TRUE}, the real shape of OCNs is plotted. If flow crosses a boundary, 
the pixel that is not contiguous to its outlet is flipped.}

\item{DrawNodes}{If \code{FALSE}, the theme is directly displayed on the river network. 
If \code{TRUE}, the theme is displayed via markers at the locations of the nodes at the RN or AG level (depending on the length of \code{theme}).
In this case, \code{NaNColor} can be used to define the color of the river network.}

\item{Node_cex}{Only effective if \code{DrawNodes == TRUE}. It sets the dimension of the markers (equivalent to 
parameter \code{cex} of function \code{\link{points}}). It can be a scalar or a vector of length \code{length(theme)}.}

\item{Node_pch}{Only effective if \code{DrawNodes == TRUE}. It sets the type of the markers (equivalent to 
parameter \code{pch} of function \code{\link{points}}). It can be a scalar or a vector of length \code{length(theme)}.}

\item{NaNColor}{Color attributed to reaches or RN nodes whose theme value is \code{NA} or \code{NaN}.}

\item{BackgroundColor}{Color used in the background of the figure. It can be either a single value, or a vector with number of components 
equal to \code{OCN$N_outlets}. If \code{OCN$N_outlets > 1} and \code{length(BackgroundColor) == OCN$N_outlets}, each color is used to identify a different catchment
(corresponding to the respective outlet). If instead \code{OCN$N_outlets > 1} and \code{length(BackgroundColor) == 1}, all catchment are colored with the same 
\code{BackgroundColor}.}

\item{add_colorbar}{If \code{TRUE}, add colorbar to the plot.}
}
\value{
No output is returned.
}
\description{
Function that draws OCNs with color of reaches or RN nodes depending on an arbitrary theme.
}
\details{
This function can be used to show how a certain spatial field varies along the river network.
}
\examples{
OCN <- aggregate_OCN(landscape_OCN(OCN_20))
draw_thematic_OCN(OCN$AG$StreamOrder, OCN, ColPalette = rainbow,
DiscreteLevels = TRUE, BackgroundColor = "#999999")
}
