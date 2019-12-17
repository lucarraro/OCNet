\name{create_OCN}
\alias{create_OCN}
\title{Create an Optimal Channel Network}
\usage{
create_OCN(dimX, dimY, nOutlet = 1, outletSide = "S",
  outletPos = round(dimX/3), periodicBoundaries = FALSE,
  typeInitialState = NULL, flowDirStart = NULL, expEnergy = 0.5,
  cellsize = 1, xllcorner = 0.5 * cellsize, yllcorner = 0.5 *
  cellsize, nIter = 30 * dimX * dimY, nUpdates = 50,
  initialNoCoolingPhase = 0, coolingRate = 50,
  showIntermediatePlots = FALSE, thrADraw = 0.002 * dimX * dimY *
  cellsize^2, saveEnergy = FALSE, saveExitFlag = FALSE,
  saveN8 = FALSE, saveN4 = FALSE, displayUpdates = 1)
}
\arguments{
\item{dimX}{Longitudinal dimension of the lattice (in number of pixels).}

\item{dimY}{Latitudinal dimension of the lattice (in number of pixels).}

\item{nOutlet}{Number of outlets. If \code{nOutlet = "All"}, all border pixels are set as outlets.}

\item{outletSide}{Side of the lattice where the outlet(s) is/are placed. 
It is a vector of characters, whose allowed values are \code{"N"} (northern side), \code{"E"}, \code{"S"}, \code{"W"}.
Its length must be equal to \code{nOutlet}.}

\item{outletPos}{Vector of positions of outlets within the sides specified by \code{outletSide}. 
If \code{outletSide[i] = "N"} or \code{"S"}, then \code{outletPos[i]} must be a natural number in the interval \code{1:dimX}; 
if \code{outletSide[i] = "W"} or \code{"E"}, then \code{outletPos[i]} must be a natural number in the interval \code{1:dimY}. 
If \code{nOutlet > 1} is specified by the user and \code{outletSide}, \code{outletPos} are not, a number of outlets equal to 
\code{nOutlet} is randomly drawn among the border pixels. Its length must be equal to \code{nOutlet}.}

\item{periodicBoundaries}{If \code{TRUE}, periodic boundaries are applied. In this case, the lattice is the planar equivalent of a torus.}

\item{typeInitialState}{Configuration of the initial state of the network. Possible values: \code{"I"} (representsting a valley); 
\code{"T"} (T-shaped drainage pattern); \code{"V"} (V-shaped drainage pattern); \code{"H"} (hip roof). Default value is set to \code{"I"}, unless when \code{nOutlet = "All"}, where default is \code{"H"}. 
See Details for explanation on initial network state in the multiple outlet case.}

\item{flowDirStart}{Matrix (dimY by dimX) with custom initial flow directions. Possible entries to \code{flowDirStart} are natural numbers between 1 and 8. 
Key is as follows (capital letters indicate cardinal directions)
\describe{
\item{\code{1}}{E (+1 column)}
\item{\code{2}}{SE (-1 row, +1 column)}
\item{\code{3}}{S (-1 row)}
\item{\code{4}}{SW (-1 row, -1 column)}
\item{\code{5}}{W (-1 column)}
\item{\code{6}}{NW (+1 row, -1 column)}
\item{\code{7}}{N (+1 row)}
\item{\code{8}}{ NE (+1 row, +1 column)}
}}

\item{expEnergy}{Exponent of the functional \code{sum(A^expEnergy)} that is minimized during the OCN search algorithm.}

\item{cellsize}{Size of a pixel in planar units.}

\item{xllcorner}{Longitudinal coordinate of the lower-left pixel.}

\item{yllcorner}{Latitudinal coordinate of the lower-left pixel.}

\item{nIter}{Number of iterations for the OCN search algorithm.}

\item{nUpdates}{Number of updates given during the OCN search process (only effective if \code{any(displayUpdates,showIntermediatePlots)=TRUE}.).}

\item{initialNoCoolingPhase, coolingRate}{Parameters of the function used to describe the temperature of the simulated annealing algorithm. See details.}

\item{showIntermediatePlots}{If \code{TRUE}, the OCN plot is updated \code{nUpdates} times during the OCN search process.
Note that, for large lattices, \code{showIntermediatePlots = TRUE} might slow down the search process considerably.}

\item{thrADraw}{Threshold drainage area value used to display the network (only effective when \code{showIntermediatePlots = TRUE}).}

\item{saveEnergy}{If \code{TRUE}, \code{Energy} is saved (see Value for its definition).}

\item{saveExitFlag}{If \code{TRUE}, \code{ExitFlag} is saved (see Value for its definition).}

\item{saveN8}{If \code{TRUE}, the adjacency matrix relative to 8-nearest-neighbours connectivity is saved.}

\item{saveN4}{If \code{TRUE}, the adjacency matrix relative to 4-nearest-neighbours connectivity is saved.}

\item{displayUpdates}{State if updates are printed on the console while the OCN search algorithm runs.
\describe{
\item{\code{0}}{No update is given.}
\item{\code{1}}{An estimate of duration is given (only if \code{dimX*dimY > 1000}, otherwise no update is given).}
\item{\code{2}}{Progress updates are given. The number of these is controlled by \code{nUpdates}}
}}
}
\value{
A list whose objects are listed below. Variables that define the network at the FD level are wrapped in the sublist \code{FD}.
Adjacency matrices describing 4- or 8- nearest-neighbours connectivity among pixels are contained in lists N4 and N8, respectively. 
\item{\code{FD$A}}{Vector (of length \code{dimX*dimY}) containing drainage area values for all FD pixels (in square planar units).}
\item{\code{FD$W}}{Adjacency matrix (\code{dimX*dimY} by \code{dimX*dimY}) at the FD level. 
It is a \code{\link{spam}} object.}
\item{\code{FD$downNode}}{Vector (of length \code{dimX*dimY}) representing the adjacency matrix at FD level in a vector form: 
if \code{FD$downNode[i] = j} then \code{FD$W[i,j] = 1}. If \code{o} is the outlet pixel, 
then \code{FD$downNode[o] = 0}.}
\item{\code{FD$X (FD$Y)}}{Vector (of length \code{dimX*dimY}) containing X (Y) coordinate values for all FD pixels.}
\item{\code{FD$nNodes}}{Number of nodes at FD level (equal to \code{dimX*dimY}).}
\item{\code{FD$outlet}}{Vector (of length \code{nOutlet}) indices of pixels at FD level corresponding to outlets.}
\item{\code{energy}}{Vector (of length \code{nIter}) of energy values for each stage of the OCN during the search algorithm (only present if \code{saveEnergy = TRUE}).}
\item{\code{exitFlag}}{Vector (of length \code{nIter}) showing the outcome of the rewiring process (only present if \code{saveExitFlag = TRUE}). Its entries can assume one of the following values:
\describe{
\item{\code{0}}{Rewiring is accepted.}
\item{\code{1}}{Rewiring is not accepted (because it does not lower \code{Energy} or according to the acceptance probability of the simulated annealing algorithm).}
\item{\code{2}}{Rewiring is invalid  because a loop in the graph was generated, therefore the network is no longer a direct acyclic graph.}
\item{\code{3}}{Rewiring is invalid because of cross-flow. This means that, for example, in a 2x2 cluster of pixel, the southwestern (SW) corner drains into the NE one, and SE drains into NW.
Although this circumstance does not imply the presence of a loop in the graph, it has no physical meaning and is thereby forbidden.}
}}
\item{\code{N4$W}}{Adjacency matrix (\code{dimX*dimY} by \code{dimX*dimY}) that describes 4-nearest-neighbours connectivity 
between pixels: \code{N4$W[i,j] = 1} if pixel \code{j} shares an edge with \code{i}, and is null otherwise. 
It is saved only if \code{saveN4 = TRUE}.} 
\item{\code{N8$W}}{Adjacency matrix (\code{dimX*dimY} by \code{dimX*dimY}) that describes 8-nearest-neighbours connectivity 
between pixels: \code{N8$W[i,j] = 1} if pixel \code{j} shares an edge or a vertex with \code{i}, and is null otherwise. 
It is saved only if \code{saveN8 = TRUE}.}  
Finally, \code{dimX}, \code{dimY}, \code{cellsize}, \code{nOutlet}, \code{periodicBoundaries}, \code{expEnergy}, 
\code{coolingRate}, \code{typeInitialState}, \code{nIter} are passed to the list as they were included in the input (except \code{nOutlet = "All"} which is converted to
 \code{2*(dimX + dimY - 2))}.
}
\description{
Function that performs the OCN search algorithm on a rectangular lattice and creates OCN at the flow direction (FD) level.
}
\details{
\emph{Simulated annealing temperature.} 
The function that expresses the temperature of the simulated annealing process is as follows: 
\describe{
 \item{if \code{i <= initialNoCoolingPhase*nIter}:}{\code{Temperature[i] = Energy[1]}}
\item{if \code{initialNoCoolingPhase*nIter < i <= nIter}:}{\code{Temperature[i] = Energy[1]*(-coolingRate*(i - InitialNocoolingPhase*nIter)/nIter)}}
 }
 where \code{i} is the index of the current iteration and \code{Energy[1] = sum(A^expEnergy)}, with \code{A} denoting 
 the vector of drainage areas corresponding to the initial state of the network. According to the simulated annealing 
 principle, a new network configuration obtained at iteration \code{i} is accepted with probability equal to 
 \code{exp((Energy[i] - Energy[i-1])/Temperature[i])} if \code{Energy[i] < Energy[i-1]}. 
 To ensure convergence, it is recommended to use \code{coolingRate} values between 20 and 200 and \code{initialNoCoolingPhase <= 0.3}. 
 Low \code{coolingRate} and high \code{initialNoCoolingPhase} values cause the network configuration to depart more significantly from the initial state.
 If \code{coolingRate < 20} is used, it is suggested to increase \code{nIter} with respect to the default value in order to guarantee convergence.

\emph{Initial network state.} 
If \code{nOutlet > 1}, the initial state is applied with regards to the outlet located at \code{outletSide[1]}, \code{outletPos[1]}. 
Subsequently, for each of the other outlets, the drainage pattern is altered within a region of maximum size \code{0.3*dimX} by \code{0.15*dimY} for outlets located at the eastern and western borders of the lattice, 
and \code{0.15*dimX} by \code{0.3*dimY} for outlets located at the southern and northern borders of the lattice. The midpoint of the long size of the regions coincides with the outlet at stake.
Within these regions, an \code{"I"}-type drainage pattern is produced if \code{typeInitialState = "I"} or \code{"T"}; a \code{"V"}-type drainage pattern is produced if \code{typeInitialState = "V"}; 
no action is performed if \code{typeInitialState = "H"}. Note that \code{typeInitialState = "H"} is the recommended choice only for large \code{nOutlet}.

\emph{Suggestions for creating "fancy" OCNs.} 
In order to generate networks spanning a realistic, non-rectangular catchment domain (in the "real-shape" view provided by \code{\link{draw_contour_OCN}}), it is convenient 
to use the option \code{periodicBoundaries = TRUE} and impose at least a couple of diagonally adjacent outlets on two opposite sides, for example \code{nOutlet = 2}, \code{outletSide = c("S", "N")}, \code{outletPos = c(1, 2)}. 
See also \code{\link{OCN_300_4out_PB}}. Note that, because the OCN search algorithm is a stochastic process, the successful generation of a "fancy" OCN is not guaranteed: indeed, it is possible that the final outcome is a 
network where most (if not all) pixels drain towards one of the two outlets, and hence such outlet is surrounded (in the "real-shape" view) by the pixels that it drains. 

In order to create a network spanning a "pear-shaped" catchment (namely where the width of the area spanned in the direction orthogonal to the main stem diminishes downstream, until it coincides with the river width at the outlet),
it is convenient to use the option \code{nOutlet = "All"} (here the value of \code{periodicBoundaries} is irrelevant) and then pick a single catchment (presumably one with rather large catchment area, see value \code{OCN$CM$A} generated by \code{\link{landscape_OCN}}) among the many generated. Note that it is not possible to predict the area spanned by such catchment \emph{a priori}. To obtain a catchment whose size is rather large compared to the size of the lattice where the 
OCN was generated, it is convenient to set \code{typeInitialState = "I"} and then pick the catchment with largest area (\code{\link{landscape_OCN}} must be run).
  
The default temperature schedule for the simulated annealing process is generally adequate for generating an OCN that does not resemble the initial network state if the size of the lattice is not too large (say, until \code{dimX*dimY <= 40000}). When \code{dimX*dimY > 40000}, it might be convenient to make use of a "warmer" temperature schedule (for example, by setting \code{coolingRate = 20} and \code{initialNoCoolingPhase = 0.15}; see also the package vignette) and/or increase \code{nIter} with respect to its default value. Note that these suggestions only pertain to the aesthetics of the final OCN; the default temperature schedule and \code{nIter} are calibrated to ensure convergence of the OCN (i.e. achievement of a local minimum of \code{Energy}, save for a reasonable threshold) also for lattices larger than \code{dimX*dimY = 40000}.
}
\examples{
# 1) creates and displays a single outlet 20x20 OCN with default options
set.seed(1)
OCN_a <- create_OCN(20, 20)
draw_simple_OCN(OCN_a)

\donttest{
# 2) creates and displays a 2-outlet OCNs with manually set outlet location, 
# and a 4-outlet OCNs with random outlet position.
set.seed(1)
old.par <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
OCN_b1 <- create_OCN(30, 30, nOutlet = 2, outletSide = c("N", "W"), outletPos = c(15, 12))
OCN_b2 <- create_OCN(30, 30, nOutlet = 4)
draw_simple_OCN(OCN_b1)
title("2-outlet OCN")
draw_simple_OCN(OCN_b2)
title("4-outlet OCN")
par(old.par)
}

\donttest{
# 3) generate 3 single-outlet OCNs on the same (100x100) domain starting from different 
# initial states, and show 20 intermediate plots and console updates.
set.seed(1)
OCN_V <- create_OCN(100, 100, typeInitialState = "V", showIntermediatePlots = TRUE, 
		nUpdates = 20, displayUpdates = 2)
OCN_T <- create_OCN(100, 100, typeInitialState = "T", showIntermediatePlots = TRUE, 
		nUpdates = 20, displayUpdates = 2)
OCN_I <- create_OCN(100, 100, typeInitialState = "I", showIntermediatePlots = TRUE, 
		nUpdates = 20, displayUpdates = 2)
}

\donttest{
# 4) generate a 2-outlet OCN and show intermediate plots. Note that different colors are used 
# to identify the two networks  (all pixels are colored because thrADraw = 0).
set.seed(1)
OCN <- create_OCN(150, 70, nOutlet = 2, outletPos = c(1, 150), outletSide = c("S", "N"),
		typeInitialState = "V", periodicBoundaries = TRUE, 
		showIntermediatePlots = TRUE, thrADraw = 0)
# The resulting networks have an irregular contour, and their outlets are located on the contour:
draw_contour_OCN(landscape_OCN(OCN))
}
}

