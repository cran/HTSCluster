\name{plot.HTSCluster}
\alias{plot.HTSCluster}
\title{Visualize results from clustering using a Poisson mixture model}
\description{
A function to visualize the clustering results obtained from a Poisson mixture model.
}
\usage{
\method{plot}{HTSCluster}(x, cluster.choice = "ICL", file.name = FALSE, 
	graphs = c("BIC.ICL", "map", "map.bycluster", "weighted.histograms", "lambda"), ...)
}
\arguments{
  \item{x}{ An object of class \code{"HTSCluster"}}
  \item{cluster.choice}{  The number of clusters desired, as selected by the "ICL"
criterion (the default), the "BIC" criterion, or an arbitrary number of clusters
among those estimated in \code{x}}
  \item{file.name}{ Optional file name if plots are to be saved in a PDF file.}
  \item{graphs}{ Type of graph to be included in plots.}
  \item{...}{Additional arguments (mainly useful for plotting)}
}
\details{
Provides the following visualizations:

1) A histogram of maximum conditional probabilities across all clusters.

2) Per-cluster boxplots of maximum conditional probabilities.

3) Weighted histograms of observation profiles (with weights equal to the 
corresponding conditional probability for each observation in each cluster),
plotted independently for each variable. Fitted densities after fitting the
Poisson mixture model are overlaid in red.

4) A global view of \eqn{\ensuremath\boldsymbol{\lambda}}{\lambda} and \eqn{\ensuremath\boldsymbol{\pi}}{\pi} values for the selected model. When the
number of conditions <= 2, bar heights represent the value of \eqn{\ensuremath\boldsymbol{\lambda}_k}{\lambda_k} for each
cluster, and bar width corresponds to the value of \eqn{\ensuremath\boldsymbol{\pi_k}}{\pi_k}. When the number of
conditions > 2, star plots are drawn for each cluster to illustrate the relative
value of \eqn{\ensuremath\boldsymbol{\lambda}_k}{\lambda_k} in each condition, where half-circle radiuses reflect the corresponding
value of \eqn{\ensuremath\boldsymbol{\pi_k}}{\pi_k}.
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

## Visualization of results (not run):

## ICL results:
## plot(run, cluster.choice = "ICL")

## BIC results:
## plot(run, cluster.choice = "BIC")

## Results for cluster size = 5:
## plot(run, cluster.choice = 5)

}
\keyword{methods}
\keyword{models}
