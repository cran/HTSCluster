\name{summary.HTSCluster}
\alias{summary.HTSCluster}
\alias{summary.HTSClusterWrapper}
\title{Summarize results from clustering using a Poisson mixture model}
\description{
A function to summarize the clustering results obtained from a Poisson mixture model.
}
\usage{
\method{summary}{HTSCluster}(object, ...)
\method{summary}{HTSClusterWrapper}(object, ...)
}
\arguments{
  \item{object}{ An object of class \code{"HTSCluster"} or \code{"HTSClusterWrapper"}}
  \item{...}{Additional arguments}
}
\details{
The summary function for an object of class \code{"HTSCluster"} provides the following summary of results:

1) Number of clusters and model selection criterion used, if applicable.

2) Number of observations across all clusters with a maximum conditional probability
greater than 90% (and corresponding percentage of total observations) for the selected
model.

3) Number of observations per cluster with a maximum conditional probability
greater than 90% (and corresponding percentage of total observations per cluster) for the
selected model.

4) \eqn{\boldsymbol{\lambda}}{\lambda} values for the selected model.

5) \eqn{\boldsymbol{\pi}}{\pi} values for the selected model.



The summary function for an object of class \code{"HTSClusterWrapper"} provides the number of clusters selected for
the BIC, ICL, DDSE, and Djump model selection approaches.
}

\references{
Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux G. (2015). Co-expression analysis of high-throughput transcriptome sequencing data with Poisson mixture models. Bioinformatics, 31(9):1420-1427.

Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011). Clustering high-throughput sequencing data with Poisson mixture models. Inria Research Report 7786. Available at \url{https://inria.hal.science/inria-00638082}.
}
\author{Andrea Rau}

\seealso{\code{\link{PoisMixClus}}, \code{\link{PoisMixClusWrapper}}}

\examples{

set.seed(12345)

## Simulate data as shown in Rau et al. (2011)
## Library size setting "A", high cluster separation
## n = 2000 observations
simulate <- PoisMixSim(n = 200, libsize = "A", separation = "high")
y <- simulate$y
conds <- simulate$conditions

## Run the PMM-II model for g = 3
## "TC" library size estimate, EM algorithm
run <- PoisMixClus(y, g = 3, 
    norm = "TC", conds = conds, init.type = "small-em")

## Summary of results:
summary(run)


}
\keyword{methods}
