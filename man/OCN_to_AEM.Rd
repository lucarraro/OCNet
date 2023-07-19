\name{OCN_to_AEM}
\alias{OCN_to_AEM}
\title{Construct asymmetric eigenvector maps (AEM) from an OCN}
\usage{
OCN_to_AEM(OCN, level = "AG", weight = NULL, resistance = "length", moranI = FALSE)
}
\arguments{
\item{OCN}{A \code{river} object as produced by \code{\link{aggregate_OCN}}.}

\item{level}{Aggregation level at which AEMs are to be calculated. It 
must be equal to either \code{"RN"} or \code{"AG"}.}

\item{weight}{Determines how and if weights should be used to compute the AEMs.
Defaults to \code{NULL}, which applies equal weights to all links of the OCN. It can
be one of \code{"gravity"}, \code{"exponential"}, \code{"linear"}, \code{"parabolic"}
or a user-specified function. See details.} 

\item{resistance}{Identifies how resisitance (i.e., the variable negatively related to 
the link weight) is calculated. Defaults to \code{"length"} (i.e., resistance equal to link length).
Alternatively, \code{resistance = "time"} uses water travel time to weigh the links.
In such a case, water velocities must be contained in the \code{OCN} (e.g. via 
\code{\link{rivergeometry_OCN}}).}

\item{moranI}{Logical. Should Moran's I statistics be computed and random tests be performed via 
\code{\link{moran.randtest}}? 
}

}
\value{
A list as produced by a call to \code{\link{aem}}. If \code{moranI = TRUE}, a \code{krandtest} resulting from
the call to \code{moran.randtest} is appended to the output list.
}
\description{
Function that computes asymmetric eigenvector maps from an OCN. These can be used as spatial variables
to assess spatial gradients in environmental or ecological data.
}
\details{
Possible character strings for \code{weight} are:
\describe{
\item{\code{"gravity"}}{\eqn{w(r) = r_{max}/r}}
\item{\code{"exponential"}}{\eqn{w(r) = \exp(-r/r_{max})}}
\item{\code{"linear"}}{\eqn{w(r) = 1 - r/r_{max}}}
\item{\code{"parabolic"}}{\eqn{w(r) = 1 - (r/r_{max})^2}}
where \eqn{w} is the weight value for a given link, \eqn{r} its resistance value and \eqn{r_{max}} the maximum resistance value across all links.
}
}
\seealso{\code{\link{aem}}, \code{\link{moran.randtest}}}

\examples{
OCN <- aggregate_OCN(landscape_OCN(OCN_20), thrA = 5)

res <- OCN_to_AEM(OCN) # unweighted AEMs
res$values # eigenvectors associates with the AEMs
plot(OCN, res$vectors[,1], drawNodes = TRUE, 
colLevels = c(-max(abs(res$vectors[,1])), max(abs(res$vectors[,1])), 100), 
	colPalette = hcl.colors(100,"Blue-Red 2")) # plot first eigenvector

res_g <- OCN_to_AEM(OCN, weight = "gravity") # weighted AEMs based on gravity model

fn <- function(r) {1 - r^0.5} 
res_f <- OCN_to_AEM(OCN, weight = fn) # weighted AEMs based on user-specified weight function

# compute Moran's I and perform permutation test to assess which eigenfunctions should be retained
\donttest{
res_g <- OCN_to_AEM(OCN, weight = "gravity", moranI = TRUE)
selectedAEM <- which(res_g$moranI$pvalue < 0.05) 
# selected eigenfunctions are those with significantly positive spatial autocorrelation

# plot selected eigenfunctions 
# (these could be e.g. used as spatial covariates in a species distribution model)
par(mfrow=c(3,4))
for (i in selectedAEM){
plot(OCN, res$vectors[,i], drawNodes = TRUE, 
	colLevels = c(-max(abs(res$vectors[,i])), max(abs(res$vectors[,i])), 100), 
	colPalette = hcl.colors(100,"Blue-Red 2"))
	title(paste0("AEM",i))
	}
}
}
