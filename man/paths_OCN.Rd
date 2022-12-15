\name{paths_OCN}
\alias{paths_OCN}
\title{Calculate paths between nodes in an Optimal Channel Network}
\usage{
paths_OCN(OCN, includePaths = FALSE, 
includeDownstreamNode = FALSE, includeUnconnectedPaths = FALSE, displayUpdates = FALSE)
}
\arguments{
\item{OCN}{List as produced by \code{\link{aggregate_OCN}}.}

\item{includePaths}{Logical. If \code{TRUE}, \code{RN$downstreamPath} and \code{AG$downstreamPath} are included to the output list. Note that this might slow down the function execution considerably, and create RAM issues for very large OCNs.}

\item{includeDownstreamNode}{Logical. If \code{TRUE}, path lengths include the length of the edge departing from the last downstream node of the path.}

\item{includeUnconnectedPaths}{Logical. If \code{TRUE}, calculate path lengths between unconnected nodes (\code{RN$downstreamLengthUnconnected} and \code{AG$downstreamLengthUnconnected}). Note that this might slow down the function execution considerably, and create RAM issues for very large OCNs.}

\item{displayUpdates}{Logical. State if updates are printed on the console while \code{landscape_OCN} runs.}
}
\value{
A list that contains all objects contained in \code{OCN}, in addition to the objects listed below. 
\item{\code{RN$downstreamPath}}{List (of length \code{OCN$RN$nNodes}) whose object \code{i} is a list (of length \code{OCN$RN$nNodes}). 
If nodes \code{i} and \code{j} are connected by a downstream path, then \code{RN$downstreamPath[[i]][[j]]} is a vector 
containing the indices of the nodes constituting such path (\code{i} and \code{j} are included). If \code{i} and \code{j} are not connected 
by a downstream path, then \code{RN$downstreamPath[[i]][[j]] = NULL}. Only present if \code{includePaths = TRUE}.}
\item{\code{RN$downstreamPathLength}}{Sparse matrix (\code{OCN$RN$nNodes} by \code{OCN$RN$nNodes}) containing length of paths between nodes 
that are connected by a downstream path; if \code{i} and \code{j} are not connected by a downstream path, then \code{RN$downstreamPathLength[i,j] = 0}. 
Note that \code{RN$downstreamPathLength[i,i] = 0} if \code{includeDownstreamNode = FALSE}; alternatively, it is \code{RN$downstreamPathLength[i,i] = OCN$RN$leng[i]}. 
It is a \code{\link{spam}} object.}
\item{\code{RN$downstreamLengthUnconnected}}{Sparse matrix (\code{OCN$RN$nNodes} by \code{OCN$RN$nNodes}). \code{RN$downstreamLengthUnconnected[i,j]} is 
the length of the downstream portion of a path joining node \code{i} to \code{j} if \code{i} and \code{j} are not connected by a downstream path.
Specifically, \code{RN$downstreamLengthUnconnected[i,j] = RN$downstreamPathLength[i,k]}, where \code{k} is a node such that there exist a downstream 
path from \code{i} to \code{k} and from \code{j} to \code{k}, and these paths are the shortest possible. Note that the length of the upstream portion 
of the path joining \code{i} to \code{j} is given by \code{RN$downstreamLengthUnconnected[j,i]}. If instead \code{i} and \code{j} are joined by a 
downstream path, then \code{RN$downstreamLengthUnconnected[i,j] = 0}. It is a \code{\link{spam}} object. Only present if \code{includeUnconnectedPaths = TRUE}.}
\item{\code{AG$downstreamPath}}{List (of length \code{OCN$AG$nNodes}) whose object \code{i} is a list (of length \code{OCN$AG$nNodes}). 
If nodes \code{i} and \code{j} are connected by a downstream path, then \code{AG$downstreamPath[[i]][[j]]} is a vector 
containing the indices of the nodes constituting such path (\code{i} and \code{j} are included). If \code{i} and \code{j} are not connected 
by a downstream path, then \code{AG$downstreamPath[[i]][[j]] = NULL}. Only present if \code{includePaths = TRUE}.}
\item{\code{AG$downstreamPathLength}}{Sparse matrix (\code{OCN$AG$nNodes} by \code{OCN$AG$nNodes}) containing length of paths between nodes 
that are connected by a downstream path; if \code{i} and \code{j} are not connected by a downstream path, then \code{AG$downstreamPathLength[i,j] = 0}. 
Note that \code{AG$downstreamPathLength[i,i] = 0} if \code{includeDownstreamNode = FALSE}; alternatively, it is \code{AG$downstreamPathLength[i,i] = OCN$AG$leng[i]}. 
It is a \code{\link{spam}} object.}
\item{\code{AG$downstreamLengthUnconnected}}{Sparse matrix (\code{OCN$AG$nNodes} by \code{OCN$AG$nNodes}). \code{AG$downstreamLengthUnconnected[i,j]} is 
the length of the downstream portion of a path joining node \code{i} to \code{j} if \code{i} and \code{j} are not connected by a downstream path.
Specifically, \code{AG$downstreamLengthUnconnected[i,j] = AG$downstreamPathLength[i,k]}, where \code{k} is a node such that there exist a downstream 
path from \code{i} to \code{k} and from \code{j} to \code{k}, and these paths are the shortest possible. Note that the length of the upstream portion 
of the path joining \code{i} to \code{j} is given by \code{AG$downstreamLengthUnconnected[j,i]}. If instead \code{i} and \code{j} are joined by a 
downstream path, then \code{AG$downstreamLengthUnconnected[i,j] = 0}. It is a \code{\link{spam}} object. Only present if \code{includeUnconnectedPaths = TRUE}.}
}
\description{
Function that determines upstream and downstream paths and path lengths between any nodes of the network at the aggregated level.
}
\examples{
# 1) Calculate paths between nodes of an OCN
OCN <- paths_OCN(aggregate_OCN(landscape_OCN(OCN_20), thrA = 4))
\dontrun{
# 2) Display distance to outlet (at the RN level) along the main stem
# of an OCN
OCN <- aggregate_OCN(landscape_OCN(OCN_250_T)) # this takes some seconds
OCN <- paths_OCN(OCN, includePaths = TRUE) # this takes some seconds

distanceToOutlet <- OCN$RN$downstreamPathLength[,OCN$RN$outlet]
farthestNode <- which(distanceToOutlet == max(distanceToOutlet))
mainStem <- OCN$RN$downstreamPath[[farthestNode]][[OCN$RN$outlet]]
theme <- rep(NaN, OCN$RN$nNodes)
theme[mainStem] <- distanceToOutlet[mainStem]

draw_thematic_OCN(theme, OCN)
title("Distance to outlet along the main stem [pixel units]")
}
}
