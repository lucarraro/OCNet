\name{paths_OCN}
\alias{paths_OCN}
\title{Calculate paths between nodes in an Optimal Channel Network}
\usage{
paths_OCN(OCN, pathsRNlevel = FALSE)
}
\arguments{
\item{OCN}{List as produced by \code{\link{aggregate_OCN}}.}

\item{pathsRNlevel}{If \code{TRUE} paths and path lengths at the RN level are copied to the output list. 
If \code{FALSE}, only path lengths at the RN level are copied to the output list.}
}
\value{
A list that contains all objects contained in \code{OCN}, in addition to the objects listed below. 
\item{\code{RN$DownstreamPath}}{(only present if \code{pathsRNlevel=TRUE}) List (of length \code{OCN$RN$Nnodes}) whose object \code{i} is a list (of length \code{OCN$RN$Nnodes}). 
If nodes \code{i} and \code{j} are connected by a downstream path, then \code{RN$DownstreamPath[[i]][[j]]} is a vector 
containing the indices of the nodes constituting such path (\code{i} and \code{j} are included). If \code{i} and \code{j} are not connected 
by a downstream path, then \code{RN$DownstreamPath[[i]][[j]] = NULL}.}
\item{\code{RN$DownstreamPathLength}}{Sparse matrix (\code{OCN$RN$Nnodes} by \code{OCN$RN$Nnodes}) containing length of paths between nodes 
that are connected by a downstream path; if \code{i} and \code{j} are not connected by a downstream path, then \code{RN$DownstreamPathLength[i,j] = 0}. 
Note that \code{RN$DownstreamPathLength[i,i] = OCN$RN$Length[i]}. It is a \code{\link{spam}} object.}
\item{\code{AG$DownstreamPath}}{List (of length \code{OCN$AG$Nnodes}) whose object \code{i} is a list (of length \code{OCN$AG$Nnodes}). 
If nodes \code{i} and \code{j} are connected by a downstream path, then \code{AG$DownstreamPath[[i]][[j]]} is a vector 
containing the indices of the nodes constituting such path (\code{i} and \code{j} are included). If \code{i} and \code{j} are not connected 
by a downstream path, then \code{AG$DownstreamPath[[i]][[j]] = NULL}.}
\item{\code{AG$DownstreamPathLength}}{Sparse matrix (\code{OCN$AG$Nnodes} by \code{OCN$AG$Nnodes}) containing length of paths between nodes 
that are connected by a downstream path; if \code{i} and \code{j} are not connected by a downstream path, then \code{AG$DownstreamPathLength[i,j] = 0}. 
Note that \code{AG$DownstreamPathLength[i,i] = OCN$AG$Length[i]}. It is a \code{\link{spam}} object.}
\item{\code{AG$DwnstrLength_unconnected}}{Sparse matrix (\code{OCN$AG$Nnodes} by \code{OCN$AG$Nnodes}). \code{AG$DwnstrLength_unconnected[i,j]} is 
the length of the downstream portion of a path joining node \code{i} to \code{j} if \code{i} and \code{j} are not connected by a downstream path.
Specifically, \code{AG$DwnstrLength_unconnected[i,j] = AG$DownstreamPathLength[i,k]}, where \code{k} is a node such that there exist a downstream 
path from \code{i} to \code{k} and from \code{j} to \code{k}, and these paths are the shortest possible. Note that the length of the upstream portion 
of the path joining \code{i} to \code{j} is given by \code{AG$DwnstrLength_unconnected[j,i]}. If instead \code{i} and \code{j} are joined by a 
downstream path, then \code{AG$DwnstrLength_unconnected[i,j] = 0}. It is a \code{\link{spam}} object.}
}
\description{
Function that determines upstream and downstream paths and path lengths between any nodes of the network at the aggregated level.
}
\examples{
# Calculate paths between nodes of an OCN
OCN <- paths_OCN(aggregate_OCN(landscape_OCN(OCN_20),A_thr=4))
}
