\name{paths_OCN}
\alias{paths_OCN}
\title{Calculate paths between nodes in an Optimal Channel Network}
\usage{
paths_OCN(OCN, pathsRN = FALSE, includeDownstreamNode = FALSE)
}
\arguments{
\item{OCN}{List as produced by \code{\link{aggregate_OCN}}.}

\item{pathsRN}{Logical. If \code{TRUE} paths and path lengths at the RN level are copied to the output list. 
If \code{FALSE}, only path lengths at the RN level are copied to the output list.}

\item{includeDownstreamNode}{Logical. If \code{TRUE}, path lengths include the length of the edge departing from the last downstream node of the path.}
}
\value{
A list that contains all objects contained in \code{OCN}, in addition to the objects listed below. 
\item{\code{RN$downstreamPath}}{(only present if \code{pathsRN = TRUE}) List (of length \code{OCN$RN$nNodes}) whose object \code{i} is a list (of length \code{OCN$RN$nNodes}). 
If nodes \code{i} and \code{j} are connected by a downstream path, then \code{RN$downstreamPath[[i]][[j]]} is a vector 
containing the indices of the nodes constituting such path (\code{i} and \code{j} are included). If \code{i} and \code{j} are not connected 
by a downstream path, then \code{RN$downstreamPath[[i]][[j]] = 0}.}
\item{\code{RN$downstreamPathLength}}{Sparse matrix (\code{OCN$RN$nNodes} by \code{OCN$RN$nNodes}) containing length of paths between nodes 
that are connected by a downstream path; if \code{i} and \code{j} are not connected by a downstream path, then \code{RN$downstreamPathLength[i,j] = 0}. 
Note that \code{RN$downstreamPathLength[i,i] = 0} if \code{includeDownstreamNode = FALSE}; alternatively, it is \code{RN$downstreamPathLength[i,i] = OCN$RN$leng[i]}. 
It is a \code{\link{spam}} object.}
\item{\code{AG$downstreamPath}}{List (of length \code{OCN$AG$nNodes}) whose object \code{i} is a list (of length \code{OCN$AG$nNodes}). 
If nodes \code{i} and \code{j} are connected by a downstream path, then \code{AG$downstreamPath[[i]][[j]]} is a vector 
containing the indices of the nodes constituting such path (\code{i} and \code{j} are included). If \code{i} and \code{j} are not connected 
by a downstream path, then \code{AG$downstreamPath[[i]][[j]] = 0}.}
\item{\code{AG$downstreamPathLength}}{Sparse matrix (\code{OCN$AG$nNodes} by \code{OCN$AG$nNodes}) containing length of paths between nodes 
that are connected by a downstream path; if \code{i} and \code{j} are not connected by a downstream path, then \code{AG$downstreamPathLength[i,j] = 0}. 
Note that \code{AG$downstreamPathLength[i,i] = 0} if \code{includeDownstreamNode = FALSE}; alternatively, it is \code{AG$downstreamPathLength[i,i] = OCN$AG$leng[i]}. 
It is a \code{\link{spam}} object.}
\item{\code{AG$downstreamLengthUnconnected}}{Sparse matrix (\code{OCN$AG$nNodes} by \code{OCN$AG$nNodes}). \code{AG$downstreamLengthUnconnected[i,j]} is 
the length of the downstream portion of a path joining node \code{i} to \code{j} if \code{i} and \code{j} are not connected by a downstream path.
Specifically, \code{AG$downstreamLengthUnconnected[i,j] = AG$downstreamPathLength[i,k]}, where \code{k} is a node such that there exist a downstream 
path from \code{i} to \code{k} and from \code{j} to \code{k}, and these paths are the shortest possible. Note that the length of the upstream portion 
of the path joining \code{i} to \code{j} is given by \code{AG$downstreamLengthUnconnected[j,i]}. If instead \code{i} and \code{j} are joined by a 
downstream path, then \code{AG$downstreamLengthUnconnected[i,j] = 0}. It is a \code{\link{spam}} object.}
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
OCN <- paths_OCN(OCN, pathsRN = TRUE) # this takes some seconds

distanceToOutlet <- OCN$RN$downstreamPathLength[,OCN$RN$outlet]
farthestNode <- which(distanceToOutlet == max(distanceToOutlet))
mainStem <- OCN$RN$downstreamPath[[farthestNode]][[OCN$RN$outlet]]
theme <- rep(NaN, OCN$RN$nNodes)
theme[mainStem] <- distanceToOutlet[mainStem]

draw_thematic_OCN(theme, OCN)
title("Distance to outlet along the main stem [pixel units]")
}
}
