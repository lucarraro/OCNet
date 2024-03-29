\name{river-class}

\docType{class}

\alias{river}
\alias{river-class}

\alias{show,river-method}
\alias{$,river-method}
\alias{$<-,river-method}
\alias{names,river-method}
\alias{[[,river,character,missing-method}
\alias{[[<-,river,character,missing-method}
\alias{[[<-,river,character,missing,ANY-method}


\title{river class}

\description{
A \code{river} object contains information on river attributes at different aggregation levels. It can represent a real river network 
(obtained via \code{rivnet::extract_river}) or an optimal channel network (obtained via \code{\link{create_OCN}}).

The content of a \code{river} object can be treated as a list, hence its objects and sublists can be accessed with both the \code{$} and \code{@} operators. 

For information on the aggregation levels and on the content of a 
\code{river} object, see \code{\link{OCNet-package}}.
}

\examples{
show(OCN_20)
names(OCN_20)

# extract or replace parts of a river object
OCN_20$dimX
OCN_20@dimX
dim <- OCN_20[["dimX"]]
OCN_20$dimX <- 1
OCN_20[["dimX"]]
OCN_20[["dimX"]] <- dim
}