\name{create_OCN}
\alias{create_OCN}
\title{Create an Optimal Channel Network}
\usage{
create_OCN(dimX, dimY, N_outlet = 1, OutletSide = "S",
  OutletPos = round(dimX/3), PeriodicBoundaries = FALSE,
  type_initialstate = NULL, FlowDirStart = NULL, ExpEnergy = 0.5,
  cellsize = 1, xllcorner = 0.5 * cellsize, yllcorner = 0.5 *
  cellsize, N_iter = 30 * dimX * dimY, N_updates = 50,
  InitialNoCoolingPhase = 0, CoolingRate = 50,
  ShowIntermediatePlots = FALSE, A_thr_draw = 0.002 * dimX * dimY *
  cellsize^2, SaveEnergy = FALSE, SaveExitFlag = FALSE,
  SaveN8 = FALSE, SaveN4 = FALSE, DisplayUpdates = 1)
}
\arguments{
\item{dimX}{Longitudinal dimension of the lattice (in number of pixels).}

\item{dimY}{Latitudinal dimension of the lattice (in number of pixels).}

\item{N_outlet}{Number of outlets. If \code{N_outlet = "All"}, all border pixels are set as outlets.}

\item{OutletSide}{Side of the lattice where the outlet(s) is/are placed. 
It is a vector of characters, whose allowed values are \code{"N"} (northern side), \code{"E"}, \code{"S"}, \code{"W"}.}

\item{OutletPos}{Vector of positions of outlets within the sides specified by \code{OutletSide}. 
If \code{OutletSide[i] = "N"} or \code{"S"}, then \code{OutletPos[i]} must be a natural number in the interval \code{1:dimX}; 
if \code{OutletSide[i] = "W"} or \code{"E"}, then \code{OutletPos[i]} must be a natural number in the interval \code{1:dimY}. 
If \code{N_outlet > 1} is specified by the user and \code{OutletSide}, \code{OutletPos} are not, a number of outlets equal to 
\code{N_outlet} is randomly drawn among the border pixels.}

\item{PeriodicBoundaries}{If \code{TRUE}, periodic boundaries are applied.}

\item{type_initialstate}{Configuration of the initial state of the network. Possible values: \code{"I"} (representsting a valley); 
\code{"T"} (T-shaped drainage pattern); \code{"V"} (V-shaped drainage pattern); \code{"H"} (hip roof). Default value is set to \code{"I"}, unless when \code{N_outlet = "All"}, where default is \code{"H"}. 
See Details for explanation on initial network state in the multiple outlet case.}

\item{FlowDirStart}{Matrix (dimY by dimX) with custom initial flow directions. Possible entries to \code{FlowDirStart} are natural numbers between 1 and 8. 
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

\item{ExpEnergy}{Exponent of the functional \code{sum(A^ExpEnergy)} that is minimized during the OCN search algorithm.}

\item{cellsize}{Size of a pixel in planar units.}

\item{xllcorner}{Longitudinal coordinate of the lower-left pixel.}

\item{yllcorner}{Latitudinal coordinate of the lower-left pixel.}

\item{N_iter}{Number of iterations for the OCN search algorithm.}

\item{N_updates}{Number of updates given during the OCN search process (only effective if \code{any(DisplayUpdates,ShowIntermediatePlots)=TRUE}.).}

\item{InitialNoCoolingPhase, CoolingRate}{Parameters of the function used to describe the temperature of the simulated annealing algorithm. See details.}

\item{ShowIntermediatePlots}{If \code{TRUE}, the OCN plot is updated \code{N_updates} times during the OCN search process.
Note that, for large lattices, \code{ShowIntermediatePlots = TRUE} might slow down the search process considerably.}

\item{A_thr_draw}{Threshold drainage area value used to display the network (only effective when \code{ShowIntermediatePlots = TRUE}).}

\item{SaveEnergy}{If \code{TRUE}, \code{Energy} is saved (see Value for its definition).}

\item{SaveExitFlag}{If \code{TRUE}, \code{ExitFlag} is saved (see Value for its definition).}

\item{SaveN8}{If \code{TRUE}, the adjacency matrix relative to 8-nearest-neighbours connectivity is saved.}

\item{SaveN4}{If \code{TRUE}, the adjacency matrix relative to 4-nearest-neighbours connectivity is saved.}

\item{DisplayUpdates}{State if updates are printed on the console while the OCN search algorithm runs.
\describe{
\item{\code{0}}{No update is given.}
\item{\code{1}}{An estimate of duration is given (only if \code{dimX*dimY > 1000}, otherwise no update is given).}
\item{\code{2}}{Progress updates are given. The number of these is controlled by \code{N_updates}}
}}
}
\value{
A list whose objects are listed below. Variables that define the network at the FD level are wrapped in the sublist \code{FD}.
Adjacency matrices describing 4- or 8- nearest-neighbours connectivity among pixels are contained in lists N4 and N8, respectively. 
\item{\code{FD$A}}{Vector (of length \code{dimX*dimY}) containing drainage area values for all FD pixels (in square planar units).}
\item{\code{FD$W}}{Adjacency matrix (\code{dimX*dimY} by \code{dimX*dimY}) at the FD level. 
It is a \code{\link{spam}} object.}
\item{\code{FD$DownNode}}{Vector (of length \code{dimX*dimY}) representing the adjacency matrix at FD level in a vector form: 
if \code{FD$DownNode[i] = j} then \code{FD$W[i,j] = 1}. If \code{o} is the outlet pixel, 
then \code{FD$DownNode[o] = 0}.}
\item{\code{FD$X (FD$Y)}}{Vector (of length \code{dimX*dimY}) containing X (Y) coordinate values for all FD pixels.}
\item{\code{FD$Nnodes}}{Number of nodes at FD level (equal to \code{dimX*dimY}).}
\item{\code{FD$Outlet}}{Vector (of length \code{N_outlet}) indices of pixels at FD level corresponding to outlets.}
\item{\code{Energy}}{Vector (of length \code{N_iter}) of energy values for each stage of the OCN during the search algorithm.}
\item{\code{ExitFlag}}{Vector (of length \code{N_iter}) showing the outcome of the rewiring process. Its entries can assume one of the following values:
\describe{
\item{\code{0}}{Rewiring is accepted.}
\item{\code{1}}{Rewiring is not accepted (because it does not lower \code{Energy} or according to the acceptance probability of the simulated annealing algorithm).}
\item{\code{2}}{Rewiring is invalid  because a loop in the graph was generated, therefore the network is no longer a direct acyclic graph.}
\item{\code{3}}{Rewiring is invalid because of cross-flow. This means that, for example, in a 2x2 cluster of pixel, the southwestern (SW) corner drains into the NE one, and SE drains into NW.
Although this circumstance does not imply the presence of a loop in the graph, it has no physical meaning and is thereby forbidden.}
}}
\item{\code{N4$W}}{Adjacency matrix (\code{dimX*dimY} by \code{dimX*dimY}) that describes 4-nearest-neighbours connectivity 
between pixels: \code{N4$W[i,j] = 1} if pixel \code{j} shares an edge with \code{i}, and is null otherwise. 
It is saved only if \code{SaveN4 = TRUE}.} 
\item{\code{N8$W}}{Adjacency matrix (\code{dimX*dimY} by \code{dimX*dimY}) that describes 8-nearest-neighbours connectivity 
between pixels: \code{N8$W[i,j] = 1} if pixel \code{j} shares an edge or a vertex with \code{i}, and is null otherwise. 
It is saved only if \code{SaveN8 = TRUE}.}  
Finally, \code{dimX}, \code{dimY}, \code{cellsize}, \code{N_outlet}, \code{PeriodicBoundaries}, \code{ExpEnergy}, 
\code{CoolingRate}, \code{type_initialstate}, \code{N_iter} are passed to the list as they were included in the input (except \code{N_outlet = "All"} which is converted to
 \code{2*(dimX + dimY - 2))}.
}
\description{
Function that performs the OCN search algorithm on a rectangular lattice and creates OCN at the flow direction (FD) level.
}
\details{
\emph{Simulated annealing temperature.} The function that expresses the temperature of the simulated annealing algorithm is as follows: 
\describe{
 \item{if \code{i <= InitialNoCoolingPhase*N_iter}:}{\code{Temperature[i] = Energy[1]}}
\item{if \code{InitialNoCoolingPhase*N_iter < i <= N_iter}:}{\code{Temperature[i] = Energy[1]*(-CoolingRate*(i - InitialNocoolingPhase*N_iter)/N_iter)}}
 }
 where \code{i} is the index of the current iteration and \code{Energy[1] = sum(A^ExpEnergy)}, with \code{A} denoting 
 the vector of drainage areas corresponding to the initial state of the network. According to the simulated annealing 
 principle, a new network configuration obtained at iteration \code{i} is accepted with probability equal to 
 \code{exp((Energy[i] - Energy[i-1])/Temperature[i])} if \code{Energy[i] < Energy[i-1]}. 
 To ensure convergence, it is recommended to use \code{CoolingRate} values between 20 and 200 and \code{InitialNoCoolingPhase <= 0.3}. 
 Low \code{CoolingRate} and high \code{InitialNoCoolingPhase} values cause the network configuration to depart more significantly from the initial state.
 If \code{CoolingRate < 20} is used, it is suggested to increase \code{N_iter} with respect to the default value in order to guarantee convergence.

\emph{Initial network state.} If \code{N_outlet > 1}, the initial state is applied with regards to the outlet located at \code{OutletSide[1]}, \code{OutletPos[1]}. 
Subsequently, for each of the other outlets, the drainage pattern is altered within a region of maximum size \code{0.3*dimX} by \code{0.15*dimY} for outlets located at the eastern and western borders of the lattice, 
and \code{0.15*dimX} by \code{0.3*dimY} for outlets located at the southern and northern borders of the lattice. The midpoint of the long size of the regions coincides with the outlet at stake.
Within these regions, an \code{"I"}-type drainage pattern is produced if \code{type_initialstate = "I"} or \code{"T"}; a \code{"V"}-type drainage pattern is produced if \code{type_initialstate = "V"}; 
no action is performed if \code{type_initialstate = "H"}. Note that \code{type_initialstate = "H"} is the recommended choice only for large \code{N_outlet}.
}
\examples{
# creates a single outlet 20x20 OCN with default options
OCN <- create_OCN(20, 20) 
}
