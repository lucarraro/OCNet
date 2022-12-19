\name{create_general_contour_OCN}
\alias{create_general_contour_OCN}
\title{Create an Optimal Channel Network on a general contour}
\usage{
create_general_contour_OCN(flowDirStart, expEnergy = 0.5,
  cellsize = 1, xllcorner = 0.5 * cellsize, yllcorner = 0.5 *
  cellsize, nIter = NULL, nUpdates = 50,
  initialNoCoolingPhase = 0, coolingRate = 1,
  showIntermediatePlots = FALSE, thrADraw = NULL, 
  easyDraw = NULL, saveEnergy = FALSE, saveExitFlag = FALSE,
  displayUpdates = 1)
}
\arguments{

\item{flowDirStart}{Matrix with custom initial flow directions. Cells outside the catchment must have value equal to \code{NaN}. 
Cells within the catchment have natural numbers between 1 and 8. Key is as follows:
\describe{
\item{\code{1}}{+1 column}
\item{\code{2}}{-1 row, +1 column}
\item{\code{3}}{-1 row}
\item{\code{4}}{-1 row, -1 column}
\item{\code{5}}{-1 column}
\item{\code{6}}{+1 row, -1 column}
\item{\code{7}}{+1 row}
\item{\code{8}}{+1 row, +1 column}
}}

\item{expEnergy}{Exponent of the functional \code{sum(A^expEnergy)} that is minimized during the OCN search algorithm.}

\item{cellsize}{Size of a pixel in planar units.}

\item{xllcorner}{Longitudinal (column-wise) coordinate of the lower-left pixel of \code{flowDirStart}.}

\item{yllcorner}{Latitudinal (row-wise) coordinate of the lower-left pixel of \code{flowDirStart}.}

\item{nIter}{Number of iterations for the OCN search algorithm. Default is 40 times the number of non-\code{NaN} pixels of \code{flowDirStart}.}

\item{nUpdates}{Number of updates given during the OCN search process (only effective if \code{any(displayUpdates,showIntermediatePlots)=TRUE}.).}

\item{initialNoCoolingPhase, coolingRate}{Parameters of the function used to describe the temperature of the simulated annealing algorithm. See details.}

\item{showIntermediatePlots}{If \code{TRUE}, the OCN plot is updated \code{nUpdates} times during the OCN search process.
Note that, for large lattices, \code{showIntermediatePlots = TRUE} might slow down the search process considerably (especially when \code{easyDraw = FALSE}).}

\item{thrADraw}{Threshold drainage area value used to display the network (only effective when \code{showIntermediatePlots = TRUE}). Default value is \code{0.002*cellsize^2*nNodes}, where \code{nNodes} is the number of non-\code{NaN} pixels of \code{flowDirStart}.}

\item{easyDraw}{Logical. If \code{TRUE}, the whole network is displayed (when \code{showIntermediatePlots = TRUE}), and pixels with drainage area lower than \code{thrADraw} are displayed in light gray. If \code{FALSE}, only pixels with drainage area greater or equal to \code{thrADraw} are displayed. Default is \code{FALSE} if \code{dimX*dimY <= 40000}, and \code{TRUE} otherwise. Note that setting \code{easyDraw = FALSE} for large networks might slow down the process considerably.}

\item{saveEnergy}{If \code{TRUE}, \code{energy} is saved (see Value for its definition).}

\item{saveExitFlag}{If \code{TRUE}, \code{exitFlag} is saved (see Value for its definition).}

\item{displayUpdates}{State if updates are printed on the console while the OCN search algorithm runs.
\describe{
\item{\code{0}}{No update is given.}
\item{\code{1}}{An estimate of duration is given (only if \code{dimX*dimY > 1000}, otherwise no update is given).}
\item{\code{2}}{Progress updates are given. The number of these is controlled by \code{nUpdates}}
}}
}
\value{
A list as in \code{\link{create_OCN}}, to which the reader is referred for detailed documentation.
However, note that in this case \code{dimX} and \code{dimY} are equal to the number of columns and rows of \code{flowDirStart}, respectively, 
while \code{nNodes} is the number of non-\code{NaN} pixels in \code{flowDirStart}. Hence, \code{nNodes} is generally lower than
\code{dimX*dimY}.
}
\description{
Function that performs the OCN search algorithm on a general contour and creates OCN at the flow direction (FD) level.
}
\examples{
OCN <- create_general_contour_OCN(flowDir)
}

