\name{plot.HTSCluster}
\alias{plot.HTSCluster}
\alias{plot.HTSClusterWrapper}
\title{Visualize results from clustering using a Poisson mixture model}
\description{
A function to visualize the clustering results obtained from a Poisson mixture model.
}
\usage{
\method{plot}{HTSCluster}(x, file.name = FALSE, 
    graphs = c("map", "map.bycluster", "lambda"), data=NA, ...)
\method{plot}{HTSClusterWrapper}(x, file.name = FALSE,
    graphs = c("capushe", "ICL", "BIC"), capushe.validation=NA, ...)
}
\arguments{
  \item{x}{ An object of class \code{"HTSCluster"} or \code{"HTSClusterWrapper"}}
  \item{file.name}{ Optional file name if plots are to be saved in a PDF file.}
  \item{graphs}{ Type of graph to be included in plots. May be equal to
  \code{"map"}, \code{"may.bycluster"}, \code{"weighted.histograms"}, and/or \code{"lambda"} for objects
  of class \code{"HTSCluster"} and \code{c("ICL", "BIC")} for objects of class
  \code{"HTSClusterWrapper"}}
  \item{capushe.validation}{ Optional number of clusters to use for capushe validation (should be less than the maximum number of clusters
specificed in the \code{"HTSClusterWrapper"} object).}
  \item{data}{
(\emph{n} x \emph{q}) matrix of observed counts for \emph{n} observations and \emph{q} variables (only required for the plotting of weighted histograms)
}
  \item{...}{Additional arguments (mainly useful for plotting)}
}
\details{
For objects of class \code{"HTSCluster"}, the plotting function provides the possibility for the following visualizations:

1) A histogram of maximum conditional probabilities across all clusters.

2) Per-cluster boxplots of maximum conditional probabilities.

3) Weighted histograms of observation profiles (with weights equal to the 
corresponding conditional probability for each observation in each cluster),
plotted independently for each variable. Fitted densities after fitting the
Poisson mixture model are overlaid in red.

4) A global view of \eqn{\ensuremath\boldsymbol{\lambda}}{\lambda} and \eqn{\ensuremath\boldsymbol{\pi}}{\pi} values for the selected model. When the
number of conditions <= 2, bar heights represent the value of \eqn{\ensuremath\boldsymbol{\lambda}_k}{\lambda_k} for each
cluster, and bar width corresponds to the value of \eqn{\ensuremath\boldsymbol{\pi_k}}{\pi_k}. 


For objects of class \code{"HTSClusterWrapper"}, the plotting function provides the possibility for one or all of the following visualizations:

1) ICL plot for all fitted models.

2) BIC plot for all fitted models.

5) Capushe diagnostic plots.
}

\references{

Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux G. (2015). Co-expression analysis of high-throughput transcriptome sequencing data with Poisson mixture models. Bioinformatics, 31(9):1420-1427.

Andrea Rau, Gilles Celeux, Marie-Laure Martin-Magniette, and Cathy Maugis-Rabusseau (2011). 
Clustering high-throughput sequencing data with Poisson mixture models. \emph{Technical report}
RR-7786, Inria Saclay -- Ile-de-France.
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

## Visualization of results (not run):
## plot(run)


}
\keyword{methods}
\keyword{models}
