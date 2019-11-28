\name{aggregate_OCN}
\alias{aggregate_OCN}
\title{Aggregate an Optimal Channel Network}
\usage{
aggregate_OCN(OCN, A_thr = 0.002 * OCN$dimX * OCN$dimY *
  OCN$cellsize^2, StreamOrderType = "Strahler", MaxReachLength = Inf)
}
\arguments{
\item{OCN}{List as produced by \code{\link{landscape_OCN}}.}

\item{A_thr}{Threshold value on drainage area used to derive the aggregated network. If \code{A_thr = 0}, no aggregation is performed: 
every FD node is also a node at the RN and AG levels. In this case, the function \code{\link{aggregate_OCN}} can still be used to compute statistics 
such as \code{OCN$AG$Upstream}, which are necessary to run the subsequent function \code{\link{paths_OCN}}.}

\item{StreamOrderType}{If \code{"Strahler"}, Strahler stream order is computed; if \code{"Shreve"}, Shreve stream order is computed.}

\item{MaxReachLength}{Maximum reach length allowed (in planar units). If the path length between a channel head and the downstream confluence 
is higher than \code{MaxReachLength}, the reach starting from the channel head will have a length up to \code{MaxReachLength}, while the next downstream 
pixel is considered as a new channel head, from which a new reach departs.}
}
\value{
A list that contains all objects contained in \code{OCN}, in addition to the objects listed below. 
New sublists RN, AG, SC, containing variables at the corresponding aggregation levels, are created.
Refer to section 4.2 of the \href{../doc/OCNet.html}{vignette} for a more detailed explanation on values \code{OCN$XX$to_YY}, where \code{XX} and \code{YY} are two random aggregation levels.
\item{\code{FD$to_RN}}{Vector (of length \code{OCN$FD$Nnodes}) whose values are equal to 0 if the FD pixel does not belong to the RN level. 
If \code{OCN$FD$to_RN[i] != 0}, then \code{OCN$FD$to_RN[i]} is the index at the RN level of the pixel whose index at the FD level is \code{i}. 
\code{OCN$FD$to_RN[i] = j} implies \code{OCN$RN$to_FD[j] = i}.}
\item{\code{FD$to_SC}}{Vector (of length \code{OCN$FD$Nnodes}) of SC indices for all nodes at the FD level. 
If \code{OCN$FD$to_SC[i] = j}, then \code{i \%in\% OCN$SC$to_FD[[j]] = TRUE}.}
\item{\code{RN$A}}{Vector (of length \code{RN$Nnodes}) containing drainage area values for all RN pixels (in square planar units).}
\item{\code{RN$W}}{Adjacency matrix (\code{RN$Nnodes} by \code{RN$Nnodes}) at the RN level. It is a \code{\link{spam}} 
object.}
\item{\code{RN$DownNode}}{Vector (of length \code{RN$Nnodes}) representing the adjacency matrix at RN level in a vector form: 
if \code{RN$DownNode[i] = j} then \code{RN$W[i,j] = 1}. If \code{o} is the outlet pixel, then \code{RN$DownNode[o] = 0}.}
\item{\code{RN$DrainageDensity}}{Drainage density of the river network, calculated as total length of the river network divided by area of the lattice. It is expressed in planar units\eqn{^{-1}}.} 
\item{\code{RN$Length}}{Vector (of length \code{RN$Nnodes}) of pixel lengths at RN level.}
\item{\code{RN$Nnodes}}{Number of nodes at the RN level.}
\item{\code{RN$Nupstream}}{Vector (of length \code{RN$Nnodes}) providing the number of nodes upstream of each node (the node itself is included).}
\item{\code{RN$Outlet}}{Vector (of length \code{OCN$FD$N_outlet}) indices of pixels at RN level corresponding to outlets.}
\item{\code{RN$Slope}}{Vector (of length \code{RN$Nnodes}) of pixel slopes at RN level.}
\item{\code{RN$to_FD}}{Vector (of length \code{RN$Nnodes}) with indices at FD level of pixels belonging to RN level. \code{RN$to_FD[i] = j} implies \code{OCN$FD$to_RN[j] = i}.}
\item{\code{RN$to_AG}}{Vector (of length \code{RN$Nnodes}) with indices at AG level of pixels belonging to RN level.}
\item{\code{RN$to_CM}}{Vector (of length \code{RN$Nnodes}) with catchment index values for each RN node. Example: \code{RN$to_CM[i] = j} if node
\code{i} drains into the outlet whose location is defined by \code{OutletSide[j]}, \code{OutletPos[j]}.}
\item{\code{RN$Upstream}}{List (of length \code{RN$Nnodes}) whose object \code{i} is a vector (of length \code{RN$Nupstream[i]}) containing 
the indices of nodes upstream of a node \code{i} (including \code{i}).}
\item{\code{RN$X}}{Vector (of length \code{RN$Nnodes}) of X coordinate of pixels at RN level.}
\item{\code{RN$Y}}{Vector (of length \code{RN$Nnodes}) of Y coordinate of pixels at RN level.}
\item{\code{RN$Z}}{Vector (of length \code{RN$Nnodes}) of Z coordinate of pixels at RN level.}
\item{\code{AG$A}}{Vector (of length \code{AG$Nnodes}) containing drainage area values for all nodes at AG level. 
If \code{i} is a channel head, then \code{AG$A[RN$to_AG[i]] = RN$A[i]}.} 
\item{\code{AG$Areach}}{Vector (of length \code{AG$Nnodes}) containing drainage area values computed by accounting for the areas drained by edges departing from AG nodes. 
In other words, \code{AG$Areach[i]} is equal to the drainage area of the last downstream node belonging to the reach that departs from \code{i} 
(namely \code{AG$Areach[i] = max(A_RN[RN$to_AG == i])}).}
\item{\code{AG$W}}{Adjacency matrix (\code{AG$Nnodes} by \code{AG$Nnodes}) at the AG level. It is a \code{\link{spam}} object.}
\item{\code{AG$DownNode}}{Vector (of length \code{AG$Nnodes}) representing the adjacency matrix at AG level in a vector form: 
if \code{AG$DownNode[i] = j} then \code{AG$W[i,j] = 1}. If \code{o} is the outlet pixel, then \code{AG$DownNode[o] = 0}.}
\item{\code{AG$Length}}{Vector (of length \code{AG$Nnodes}) of lengths of edges departing from nodes at AG level. 
If \code{o} is an outlet node (i.e. (\code{o \%in\% AG$Outlet) = TRUE}), then \code{AG$Length[i] = 0}.}
\item{\code{AG$Nnodes}}{Number of nodes resulting from the aggregation process.}
\item{\code{AG$Nupstream}}{Vector (of length \code{AG$Nnodes}) providing the number of nodes upstream of each node (the node itself is included).}
\item{\code{AG$Outlet}}{Vector (of length \code{OCN$FD$N_outlet}) with indices of outlet nodes, i.e. nodes whose \code{AG$DownNode} value is 0.} 
\item{\code{AG$Slope}}{Vector (of length \code{AG$Nnodes}) of slopes at AG level. It represents the (weighted) average slope of edges departing from nodes. 
If \code{i} is an outlet node (i.e. (\code{i \%in\% AG$Outlet) = TRUE}), then \code{AG$Slope[i] = NaN}.}
\item{\code{AG$StreamOrder}}{Vector (of length \code{AG$Nnodes}) of stream order values for each node. If \code{StreamOrderType = "Strahler"}, 
Strahler stream order is computed. If \code{StreamOrderType = "Shreve"}, Shreve stream order is computed.}
\item{\code{AG$Upstream}}{List (of length \code{AG$Nnodes}) whose object \code{i} is a vector (of length \code{AG$Nupstream[i]}) containing 
the indices of nodes upstream of a node \code{i} (including \code{i}).}
\item{\code{AG$to_FD}}{List (of length \code{AG$Nnodes}) whose object \code{i} is a vector of indices of FD pixels constituting the edge 
departing from node \code{i}. If \code{i} is an outlet node, then \code{AG$to_FD} is empty.}
\item{\code{AG$to_RN}}{List (of length \code{AG$Nnodes}) whose object \code{i} is a vector of indices of RN pixels constituting the edge 
departing from node \code{i}. If \code{i} is an outlet node, then \code{AG$to_RN} is empty.}
\item{\code{AG$to_CM}}{Vector (of length \code{AG$Nnodes}) with catchment index values for each AG node. Example: \code{AG$to_CM[i] = j} if node
\code{i} drains into the outlet whose location is defined by \code{OutletSide[j]}, \code{OutletPos[j]}.}
\item{\code{AG$X}}{Vector (of length \code{AG$Nnodes}) of X coordinates (in planar units) of nodes at the AG level. 
These correspond to the X coordinates of the pixels constituting the channel heads. If \code{i} is an outlet node, then \code{AG$X[i] = NaN}.}
\item{\code{AG$Y}}{Vector (of length \code{AG$Nnodes}) of Y coordinates (in planar units) of nodes at the AG level. 
These correspond to the Y coordinates of the pixels constituting the channel heads. If \code{i} is an outlet node, then \code{AG$Y[i] = NaN}.}
\item{\code{AG$Z}}{Vector (of length \code{AG$Nnodes}) of elevation values (in elevational units) of nodes at the AG level. 
These correspond to the elevations of the pixels constituting the channel heads. If \code{i} is an outlet node, then \code{AG$Z[i] = NaN}.}
\item{\code{SC$Alocal}}{Vector (of length \code{SC$Nnodes}) with values of subcatchment area, that is the number of FD pixels 
(multiplied by \code{OCN$FD$cellsize^2}) that constitutes a subcatchment. If \code{o} is an outlet node, then \code{Alocal[o] = 0}.}
\item{\code{SC$W}}{Adjacency matrix (\code{SC$Nnodes} by \code{SC$Nnodes}) at the subcatchment level. Two subcatchments are connected if they share a border. 
Note that this is not a flow connection. Unlike the adjacency matrices at levels FD, RN, AG, this matrix is symmetric. It is a \code{\link{spam}} object.
If \code{o} is an outlet node, then \code{SC$W[o,]} and \code{SC$W[,o]} only contain zeros (i.e., \code{o} is unconnected to the other nodes).}
\item{\code{SC$Nnodes}}{Number of subcatchments into which the lattice is partitioned. If \code{N_outlet = 1}, then \code{SC$Nnodes = AG$Nnodes}. 
If multiple outlets are present, \code{SC$Nnodes} might be greater than \code{AG$Nnodes} in the case when some catchments have drainage area lower than \code{A_thr}. 
In this case, the indices from \code{AG$Nnodes + 1} to \code{SC$Nnodes} identify subcatchment that do not have a corresponding AG node.}
\item{\code{SC$to_FD}}{List (of length \code{SC$Nnodes}) whose object \code{i} is a vector of indices of FD pixels constituting the subcatchment \code{i}.}
\item{\code{SC$X}}{Vector (of length \code{SC$Nnodes}) of X coordinates of subcatchment centroids. If \code{o} is an outlet node, then \code{SC$X[o] = NaN}.}
\item{\code{SC$Y}}{vector (of length \code{SC$Nnodes}) of Y coordinates of subcatchment centroids. If \code{o} is an outlet node, then \code{SC$Y[o] = NaN}.}
\item{\code{SC$Z}}{vector (of length \code{SC$Nnodes}) of average subcatchment elevation. If \code{o} is an outlet node, then \code{SC$Z[o] = NaN}.}
Finally, \code{A_thr} is added to the list.
}
\description{
Function that, given an OCN, builds the network at the river network (RN), aggregated (AG), subcatchment (SC), and catchment (CM) levels.
}
\details{
Note that each node (and the corresponding edge exiting from it, in the case of non-outlet nodes) at the AG level corresponds to 
a subcatchment at the SC level that shares the same index: for instance, \code{SC$to_FD[i]} contains all elements of 
\code{AG$to_FD[i]} (that is, the indices of pixels at FD level that constitute the edge departing from node \code{i} 
are also part of subcatchment \code{i}).
}
\examples{
OCN <- aggregate_OCN(landscape_OCN(OCN_20))
}
