\name{summary.HTSCluster}
\alias{summary.HTSCluster}
\title{Summarize results from clustering using a Poisson mixture model}
\description{
A function to summarize the clustering results obtained from a Poisson mixture model.
}
\usage{
\method{summary}{HTSCluster}(object, cluster.choice = "ICL", ...)
}
\arguments{
  \item{object}{ An object of class \code{"HTSCluster"}}
  \item{cluster.choice}{ The number of clusters desired, as selected by the "ICL"
criterion (the default), the "BIC" criterion, or an arbitrary number of clusters
among those estimated in \code{x}}
  \item{...}{Additional arguments}
}
\details{
Provides the following summary of results:

1) Selected number of clusters (for all \code{cluster.choice}) and BIC or ICL
values for \code{cluster.choice = "BIC"} and \code{cluster.choice = "ICL"}, respectively.

2) Number of observations across all clusters with a maximum conditional probability
greater than 90% (and corresponding percentage of total observations) for the selected
model.

3) Number of observations per cluster with a maximum conditional probability
greater than 90% (and corresponding percentage of total observations per cluster) for the
selected model.

4) \eqn{\ensuremath\boldsymbol{\lambda}}{\lambda} values for the selected model.

5) \eqn{\ensuremath\boldsymbol{\pi}}{\pi} values for the selected model.

}

\references{
Andrea Rau, Gilles Celeux, Marie-Laure Martin-Magniette, and Cathy Maugis-Rabusseau (2011). 
Clustering high-throughput sequencing data with Poisson mixture models. \emph{Technical report}
RR-7786, Inria Saclay -- Ile-de-France.
}
\author{Andrea Rau}

\seealso{\code{\link{PoisMixClus}}}

\examples{

set.seed(12345)

## Simulate data as shown in Rau et al. (2011)
## Library size setting "A", high cluster separation
## n = 2000 observations
simulate <- PoisMixSim(n = 200, libsize = "A", separation = "high")
y <- simulate$y
conds <- simulate$conditions

## Run the PMM-II model for g = {3, 4, 5}
## "TC" library size estimate, EM algorithm
run <- PoisMixClus(y, gmin = 3, gmax = 5, lib.size = TRUE,
    lib.type = "TC", conds = conds, init.type = "small-em")

## Summary of results:

## ICL results:
summary(run, cluster.choice = "ICL")

## BIC results:
## summary(run, cluster.choice = "BIC")

## Results for cluster size = 5:
## summary(run, cluster.choice = 5)


}
\keyword{methods}
