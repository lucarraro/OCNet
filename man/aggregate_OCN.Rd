\name{aggregate_OCN}
\alias{aggregate_OCN}
\title{Aggregate an Optimal Channel Network}
\usage{
aggregate_OCN(OCN, thrA = 0.002 * OCN$FD$nNodes *
  OCN$cellsize^2, streamOrderType = "Strahler", maxReachLength = Inf,
  breakpoints = NULL)
}
\arguments{
\item{OCN}{List as produced by \code{\link{landscape_OCN}}.}

\item{thrA}{Threshold value on drainage area used to derive the aggregated network. If \code{thrA = 0}, no aggregation is performed: 
every FD node is also a node at the RN and AG levels. In this case, the function \code{aggregate_OCN} can still be used to compute statistics 
such as \code{OCN$AG$streamOrder}.}

\item{streamOrderType}{If \code{"Strahler"}, Strahler stream order is computed; if \code{"Shreve"}, Shreve stream order is computed.}

\item{maxReachLength}{Maximum reach length allowed (in planar units). If the path length between a channel head and the downstream confluence 
is higher than \code{maxReachLength}, the reach starting from the channel head will have a length up to \code{maxReachLength}, while the next downstream 
pixel is considered as a new channel head, from which a new reach departs.}

\item{breakpoints}{Indices of additional nodes at the RN level that should be also nodes at the AG level (beyond source, confluence, outlet nodes 
and AG nodes determined via \code{maxReachLength}). To determine such indices, a preliminary run of \code{aggregate_OCN} with the same \code{thrA} 
would be required (see example 3).}
}
\value{
A list that contains all objects contained in \code{OCN}, in addition to the objects listed below. 
New sublists RN, AG, SC, containing variables at the corresponding aggregation levels, are created.
Refer to section 4.2 of the \href{../doc/OCNet.html}{vignette} for a more detailed explanation on values \code{OCN$XX$toYY}, where \code{XX} and \code{YY} are two random aggregation levels.
\item{\code{FD$toRN}}{Vector (of length \code{OCN$FD$nNodes}) whose values are equal to 0 if the FD node is not a node at the RN level. 
If \code{FD$toRN[i] != 0}, then \code{FD$toRN[i]} is the index at the RN level of the node whose index at the FD level is \code{i}. Thereby, 
\code{FD$toRN[i] = j} implies \code{RN$toFD[j] = i}.}
\item{\code{FD$toSC}}{Vector (of length \code{OCN$FD$nNodes}) of SC indices for all nodes at the FD level. 
If \code{OCN$FD$toSC[i] = j}, then \code{i \%in\% OCN$SC$toFD[[j]] = TRUE}.}
\item{\code{RN$A}}{Vector (of length \code{RN$nNodes}) containing drainage area values for all RN nodes (in square planar units).}
\item{\code{RN$W}}{Adjacency matrix (\code{RN$nNodes} by \code{RN$nNodes}) at the RN level. It is a \code{\link{spam}} 
object.}
\item{\code{RN$downNode}}{Vector (of length \code{RN$nNodes}) representing the adjacency matrix at RN level in a vector form: 
if \code{RN$downNode[i] = j} then \code{RN$W[i,j] = 1}. If \code{o} is the outlet node, then \code{RN$downNode[o] = 0}.}
\item{\code{RN$drainageDensity}}{Drainage density of the river network, calculated as total length of the river network divided by area of the lattice. It is expressed in planar units^(-1).} 
\item{\code{RN$leng}}{Vector (of length \code{RN$nNodes}) of lengths of edges departing from nodes at the RN level. Its values are equal to either \code{0} (if the corresponding
node is an outlet), \code{OCN$cellsize} (if the corresponding flow direction is horizontal/vertical), or \code{sqrt(2)*OCN$cellsize} (diagonal flow).}
\item{\code{RN$nNodes}}{Number of nodes at the RN level.}
\item{\code{RN$nUpstream}}{Vector (of length \code{RN$nNodes}) providing the number of nodes upstream of each node (the node itself is included).}
\item{\code{RN$outlet}}{Vector (of length \code{OCN$FD$nOutlet}) indices of nodes at RN level corresponding to outlets.}
\item{\code{RN$Slope}}{Vector (of length \code{RN$nNodes}) of pixel slopes at RN level.}
\item{\code{RN$toAG}}{Vector (of length \code{RN$nNodes}) whose values are equal to 0 if the RN node is not a node at the AG level. 
If \code{RN$toAG[i] != 0}, then \code{RN$toAG[i]} is the index at the AG level of the node whose index at the RN level is \code{i}. Thereby, 
\code{RN$toAG[i] = j} implies \code{AG$toRN[j] = i}.}
\item{\code{RN$toAGReach}}{Vector (of length \code{RN$nNodes}) identifying to which edge (reach) the RN nodes belong. If \code{RN$toAGReach[i] = j}, 
the RN node \code{i} belongs to the edge departing from from the AG node \code{j} (which implies that it may correspond to the AG node \code{j} itself.)}
\item{\code{RN$toFD}}{Vector (of length \code{RN$nNodes}) with indices at FD level of nodes belonging to RN level. \code{RN$toFD[i] = j} implies \code{OCN$FD$toRN[j] = i}.}
\item{\code{RN$toCM}}{Vector (of length \code{RN$nNodes}) with catchment index values for each RN node. Example: \code{RN$toCM[i] = j} if node
\code{i} drains into the outlet whose location is defined by \code{outletSide[j]}, \code{outletPos[j]}.}
\item{\code{RN$upstream}}{List (of length \code{RN$nNodes}) whose object \code{i} is a vector (of length \code{RN$nUpstream[i]}) containing 
the indices of nodes upstream of a node \code{i} (including \code{i}).}
\item{\code{RN$X}, \code{RN$Y}}{Vectors (of length \code{RN$nNodes}) of X, Y coordinates of nodes at RN level.}
\item{\code{RN$Z}}{Vector (of length \code{RN$nNodes}) of Z coordinates of nodes at RN level.}
\item{\code{AG$A}}{Vector (of length \code{AG$nNodes}) containing drainage area values for all nodes at AG level. 
If \code{i} is a channel head, then \code{AG$A[RN$toAG[i]] = RN$A[i]}.} 
\item{\code{AG$AReach}}{Vector (of length \code{AG$nNodes}) containing drainage area values computed by accounting for the areas drained by edges departing from AG nodes. 
In other words, \code{AG$AReach[i]} is equal to the drainage area of the last downstream node belonging to the reach that departs from \code{i} 
(namely \code{AG$AReach[i] = max(RN$A[RN$toAG == i])}).}
\item{\code{AG$W}}{Adjacency matrix (\code{AG$nNodes} by \code{AG$nNodes}) at the AG level. It is a \code{\link{spam}} object.}
\item{\code{AG$downNode}}{Vector (of length \code{AG$nNodes}) representing the adjacency matrix at AG level in a vector form: 
if \code{AG$downNode[i] = j} then \code{AG$W[i,j] = 1}. If \code{o} is the outlet node, then \code{AG$downNode[o] = 0}.}
\item{\code{AG$leng}}{Vector (of length \code{AG$nNodes}) of lengths of edges departing from nodes at AG level. Note that \code{AG$leng[i] = sum(RN$leng[RN$toAG == i])}.
If \code{o} is an outlet node (i.e. (\code{o \%in\% AG$outlet) = TRUE}), then \code{AG$leng[i] = 0}.}
\item{\code{AG$nNodes}}{Number of nodes resulting from the aggregation process.}
\item{\code{AG$nUpstream}}{Vector (of length \code{AG$nNodes}) providing the number of nodes (at the AG level) upstream of each node (the node itself is included).}
\item{\code{AG$outlet}}{Vector (of length \code{OCN$FD$nOutlet}) with indices of outlet nodes, i.e. nodes whose \code{AG$downNode} value is 0.} 
\item{\code{AG$slope}}{Vector (of length \code{AG$nNodes}) of slopes at AG level. It represents the (weighted) average slope of edges departing from nodes. 
If \code{i} is an outlet node (i.e. (\code{i \%in\% AG$outlet) = TRUE}), then \code{AG$slope[i] = NaN}.}
\item{\code{AG$streamOrder}}{Vector (of length \code{AG$nNodes}) of stream order values for each node. If \code{streamOrderType = "Strahler"}, 
Strahler stream order is computed. If \code{streamOrderType = "Shreve"}, Shreve stream order is computed.}
\item{\code{AG$upstream}}{List (of length \code{AG$nNodes}) whose object \code{i} is a vector (of length \code{AG$nUpstream[i]}) containing 
the indices of nodes (at the AG level) upstream of a node \code{i} (including \code{i}).}
\item{\code{AG$toFD}}{Vector of length \code{AG$nNodes}) with with indices at FD level of nodes belonging to AG level. \code{AG$toFD[i] = j} implies \code{OCN$FD$toAG[j] = i}. }
\item{\code{AG$ReachToFD}}{List (of length \code{AG$nNodes}) whose object \code{i} is a vector of indices of FD nodes constituting the edge 
departing from node \code{i}.}
\item{\code{AG$toRN}}{Vector of length \code{AG$nNodes}) with with indices at RN level of nodes belonging to AG level. \code{AG$toRN[i] = j} implies \code{OCN$FD$toRN[j] = i}. }
\item{\code{AG$ReachToRN}}{List (of length \code{AG$nNodes}) whose object \code{i} is a vector of indices of RN nodes constituting the edge 
departing from node \code{i}.}
\item{\code{AG$toCM}}{Vector (of length \code{AG$nNodes}) with catchment index values for each AG node. Example: \code{AG$toCM[i] = j} if node
\code{i} drains into the outlet whose location is defined by \code{outletSide[j]}, \code{outletPos[j]}.}
\item{\code{AG$X}, \code{AG$Y}}{Vectors (of length \code{AG$nNodes}) of X, Y coordinates (in planar units) of nodes at the AG level. 
These correspond to the X, Y coordinates of the nodes constituting the upstream tips of the reaches. 
If \code{i} and \code{j} are such that \code{AG$X[i] == RN$X[j]} and \code{AG$Y[i] == RN$Y[j]}, then \code{AG$A[i] = RN$A[j]}.}
\item{\code{AG$XReach}, \code{AG$YReach}}{Vector (of length \code{AG$nNodes}) of X, Y coordinates (in planar units) of the downstream tips of the reaches. 
If \code{i} and \code{j} are such that \code{AG$XReach[i] == RN$X[j]} and \code{AG$YReach[i] == RN$Y[j]}, then \code{AG$AReach[i] = RN$A[j]}. 
If \code{o} is an outlet node, then \code{AG$XReach = NaN}, \code{AG$YReach = NaN}.}
\item{\code{AG$Z}}{Vector (of length \code{AG$nNodes}) of elevation values (in elevational units) of nodes at the AG level. 
These correspond to the elevations of the nodes constituting the upstream tips of the reaches.}
\item{\code{AG$ZReach}}{Vector (of length \code{AG$nNodes}) of Z coordinates (in elevational units) of the downstream tips of the reaches.
If \code{o} is an outlet node, then \code{AG$ZReach = NaN}.}
\item{\code{SC$ALocal}}{Vector (of length \code{SC$nNodes}) with values of subcatchment area, that is the number of FD pixels 
(multiplied by \code{OCN$FD$cellsize^2}) that constitutes a subcatchment. If \code{o} is an outlet node, then \code{ALocal[o] = 0}.}
\item{\code{SC$W}}{Adjacency matrix (\code{SC$nNodes} by \code{SC$nNodes}) at the subcatchment level. Two subcatchments are connected if they share a border. 
Note that this is not a flow connection. Unlike the adjacency matrices at levels FD, RN, AG, this matrix is symmetric. It is a \code{\link{spam}} object.
If \code{o} is an outlet node, then \code{SC$W[o,]} and \code{SC$W[,o]} only contain zeros (i.e., \code{o} is unconnected to the other nodes).}
\item{\code{SC$nNodes}}{Number of subcatchments into which the lattice is partitioned. If \code{nOutlet = 1}, then \code{SC$nNodes = AG$nNodes}. 
If multiple outlets are present, \code{SC$nNodes} might be greater than \code{AG$nNodes} in the case when some catchments have drainage area lower than \code{thrA}. 
In this case, the indices from \code{AG$nNodes + 1} to \code{SC$nNodes} identify subcatchment that do not have a corresponding AG node.}
\item{\code{SC$toFD}}{List (of length \code{SC$nNodes}) whose object \code{i} is a vector of indices of FD pixels constituting the subcatchment \code{i}.}
\item{\code{SC$X}, \code{SC$Y}}{Vectors (of length \code{SC$nNodes}) of X, Y coordinates (in planar units) of subcatchment centroids.}
\item{\code{SC$Z}}{Vector (of length \code{SC$nNodes}) of average subcatchment elevation (in elevational units).}
Finally, \code{thrA} is added to the list.
}
\description{
Function that, given an OCN, builds the network at the river network (RN), aggregated (AG), subcatchment (SC), and catchment (CM) levels.
}
\details{
Note that each node (and the corresponding edge exiting from it, in the case of non-outlet nodes) at the AG level corresponds to 
a subcatchment at the SC level that shares the same index: for instance, \code{SC$toFD[i]} contains all elements of 
\code{AG$toFD[i]} (that is, the indices of pixels at FD level that constitute the edge departing from node \code{i} 
are also part of subcatchment \code{i}).
}
\examples{
# 1) aggregate a 20x20 OCN by imposing thrA = 4
OCN <- aggregate_OCN(landscape_OCN(OCN_20), thrA = 4)

\donttest{
# 2) explore the effects of thrA and maxReachLength on a large OCN
OCN <- landscape_OCN(OCN_250_T) # it takes some seconds
OCN_a <- aggregate_OCN(OCN, thrA = 200) # it takes some seconds
OCN_b <- aggregate_OCN(OCN, thrA = 1000) # it takes some seconds
OCN_c <- aggregate_OCN(OCN, thrA = 1000, maxReachLength = 20) # it takes some seconds

old.par <- par(no.readonly = TRUE)
par(mfrow = c(1,3))
draw_subcatchments_OCN(OCN_a)
points(OCN_a$AG$X, OCN_a$AG$Y, pch = 19, col = "#0044bb") 
title(paste("No. AG nodes = ", as.character(OCN_a$AG$nNodes),
		sep=""))
draw_subcatchments_OCN(OCN_b)
points(OCN_b$AG$X, OCN_b$AG$Y, pch = 19, col = "#0044bb") 
title(paste("No. AG nodes = ", as.character(OCN_b$AG$nNodes),
		sep=""))
draw_subcatchments_OCN(OCN_c)
points(OCN_c$AG$X, OCN_c$AG$Y, pch = 19, col = "#0044bb") 
title(paste("No. AG nodes = ", as.character(OCN_c$AG$nNodes),
		sep=""))
par(old.par)
}

# 3) use of breakpoints
OCN <- aggregate_OCN(landscape_OCN(OCN_20), thrA = 5)
draw_thematic_OCN(OCN, drawNodes=TRUE)
# add an AG node downstream of node 1 at AG level
new_node_RN <- OCN$RN$downNode[OCN$AG$toRN[1]]
OCN2 <- aggregate_OCN(landscape_OCN(OCN_20), thrA = 5, breakpoints = new_node_RN)
draw_thematic_OCN(OCN2, drawNodes = TRUE)
points(OCN$RN$X[new_node_RN], OCN$RN$Y[new_node_RN], 
	pch = 19, col = "red") # this node has been added
}
