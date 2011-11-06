\name{PoisMixSim}
\alias{PoisMixSim}
\title{
Simulate data from a Poisson mixture model
}
\description{
This function simulates data from a Poisson mixture model, as described by Rau et al. (2011). Data are simulated with varying expression level (\eqn{w_i}) for 4 clusters. Clusters may be simulated with \dQuote{high} or \dQuote{low} separation, and three different options are available for the library size setting: \dQuote{equal}, \dQuote{A}, and \dQuote{B}, as described by Rau et al. (2011).
}
\usage{
PoisMixSim(n = 2000, libsize, separation)
}

\arguments{
  \item{n}{
Number of observations
}
  \item{libsize}{
The type of library size difference to be simulated (\dQuote{\code{equal}}, \dQuote{\code{A}}, or \dQuote{\code{B}}, as described by Rau et al. (2011))
}
  \item{separation}{
Cluster separation (\dQuote{\code{high}} or \dQuote{\code{low}}, as described by Rau et al. (2011))
}
}
\value{
\item{y }{(\emph{n} x \emph{q}) matrix of simulated counts for \emph{n} observations and \emph{q} variables}
\item{labels }{Vector of length \emph{n} defining the true cluster labels of the simulated data}
\item{pi }{Vector of length 4 (the number of clusters) containing the true value of \eqn{\ensuremath\boldsymbol{\pi}}{\pi}}
\item{lambda }{(\emph{d} x \emph{4}) matrix of \eqn{\ensuremath\boldsymbol{\lambda}}{\lambda} values for \emph{d} conditions (3 in the case of \code{libsize =} \dQuote{\code{equal}} or \dQuote{\code{A}}, and 2 otherwise) in 4 clusters (see note below)}
\item{w }{Row sums of \code{y} (estimate of \eqn{\hat{w}})}
\item{conditions }{Vector of length \emph{q} defining the condition (treatment group) for each variable (column) in \code{y}}
}
\references{
Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C (2011). Clustering high-throughput sequencing data with Poisson mixture models. Inria Research Report 7786. Available at \url{http://hal.inria.fr/inria-00638082}.
}
\author{
Andrea Rau <\url{andrea.rau@jouy.inra.fr}>
}
\note{
If one or more observations are simulated such that all variables have a value of 0, those rows are removed from the data matrix; as such, in some cases the simulated data \code{y} may have less than \code{n} rows.

The PMM-I model includes the parameter constraint \eqn{\sum_k \lambda_{jk} r_j = 1}, where \eqn{r_j} is the number of replicates in condition (treatment group) \eqn{j}. Similarly, the parameter constraint in the PMM-II model is \eqn{\sum_j \sum_l \lambda_{jk}s_{jl} = 1}, where \eqn{s_{jl}} is the library size for replicate \emph{l} of condition \emph{j}. The value of \code{lambda} corresponds to that used to generate the simulated data, where the library sizes were set as described in Table 2 of Rau et al. (2011). However, due to variability in the simulation process, the actually library sizes of the data \code{y} are not exactly equal to these values; this means that the value of \code{lambda} may not be directly compared to an estimated value of \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}} as obtained from the \code{\link{PoisMixClus}} function.
}

\examples{

set.seed(12345)

## Simulate data as shown in Rau et al. (2011)
## Library size setting "A", high cluster separation
## n = 200 observations

simulate <- PoisMixSim(n = 200, libsize = "A", separation = "high")
y <- simulate$y
conds <- simulate$conditions

}
\keyword{ datagen }

