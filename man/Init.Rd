\name{Init}
\alias{emInit}
\alias{kmeanInit}
\title{
Small EM parameter initialization
}
\description{
These functions implement the Small EM initialization strategy described in Rau et al. (2011) to obtain initial values for the parameters of a Poisson mixture model. 
}
\usage{
emInit(y, g, conds, lib.size, lib.type = "TC", alg.type = "EM", 
    starts = 5, verbose = FALSE)

kmeanInit(y, g, conds, lib.size, lib.type = "TC")
}
\arguments{
  \item{y}{
(\emph{n} x \emph{q}) matrix of observed counts for \emph{n} observations and \emph{q} variables
}
  \item{g}{
Number of clusters
}
  \item{conds}{
Vector of length \emph{q} defining the condition (treatment group) for each variable (column) in \code{y}
}
  \item{lib.size}{
If \code{FALSE}, the library size parameter is not included in the model (i.e., the PMM-I model). If \code{TRUE}, the library size parameter is included in the Poisson mixture model (i.e., the PMM-II model)
}
  \item{lib.type}{
If \code{lib.size = TRUE}, the type of estimator to be used for the library size parameter (\dQuote{\code{TC}} for total count, \dQuote{\code{Q}} for quantile, and \dQuote{\code{MedRatio}} for the median ratio of Anders and Huber (2010))
}
  \item{alg.type}{
Algorithm to be used for parameter estimation (\dQuote{\code{EM}} or \dQuote{\code{CEM}} for the EM or CEM algorithms, respectively)
}
  \item{starts}{
The number independent runs with the Small-EM algorithm (with a default value of five)
}
  \item{verbose}{
If \code{TRUE}, include verbose output
}
}
\details{
To initialize parameter values for the EM and CEM algorithms we use a so-called Small-EM strategy (Biernacki et al., 2003) using the \code{emInit} function. Five independent times, the following procedure is used to obtain parameter values: first, a K-means algorithm (MacQueen, 1967) is run to partition the data into \code{g} clusters (\eqn{\hat{\ensuremath\boldsymbol{z}}^{(0)}}{\hat{z}^(0)}). Second, initial parameter values \eqn{\ensuremath\boldsymbol{\pi}^{(0)}}{\pi^(0)} and \eqn{\ensuremath\boldsymbol{\lambda}^{(0)}}{\lambda^(0)} are calculated (see Rau et al. (2011) for details). Third, five iterations of an EM algorithm are run, using \eqn{\ensuremath\boldsymbol{\pi}^{(0)}}{\pi^(0)} and \eqn{\ensuremath\boldsymbol{\lambda}^{(0)}}{\lambda^(0)} as initial values. Finally, among the five sets of parameter values, we use \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}} and \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}} corresponding to the highest log likelihood or completed log likelihood to initialize the subsequent full EM or CEM algorithms, respectively.
}
\value{
\item{pi.init }{Vector of length \code{g} containing the estimate for \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}} corresponding to the highest log likelihood (or completed log likelihood) from the Small-EM inialization strategy. }
\item{lambda.init }{(\emph{d} x \code{g}) matrix containing the estimate of \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}} corresponding to the highest log likelihood (or completed log likelihood) from the Small-EM initialization strategy, where \emph{d} is the number of conditions and \code{g} is the number of clusters. }
}
\references{
Anders, S. and Huber, W. (2010) Differential expression analysis for sequence count data. \emph{Genome Biology}, \bold{11}(R106), 1-28.

Biernacki, C., Celeux, G., Govaert, G. (2003) Choosing starting values for the EM algorithm for getting the highest likelhiood in multivariate Gaussian mixture models. \emph{Computational Statisitcs and Data Analysis}, \bold{41}(1), 561-575.

MacQueen, J. B. (1967) Some methods for classification and analysis of multivariate observations. In \emph{Proceedings of the 5th Berkeley Symposium on Mathematical Statistics and Probability}, number 1, pages 281-297. Berkeley, University of California Press.

Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C (2011). Clustering high-throughput sequencing data with Poisson mixture models. Inria Research Report 7786. Available at \url{http://hal.inria.fr/inria-00638082}.
}
\author{
Andrea Rau <\url{andrea.rau@jouy.inra.fr}>
}

\seealso{
\code{\link{PoisMixClus} for Poisson mixture model estimation and model selection}
}
\examples{

set.seed(12345)

## Simulate data as shown in Rau et al. (2011)
## Library size setting "A", high cluster separation
## n = 500 observations

simulate <- PoisMixSim(n = 500, libsize = "A", separation = "high")
y <- simulate$y
conds <- simulate$conditions

## Calculate initial values for lambda and pi using the Small-EM
## initialization (4 classes, PMM-II model with "TC" library size)

init.values <- emInit(y, g = 4, conds, lib.size = TRUE, 
    lib.type = "TC", alg.type = "EM")
pi.init <- init.values$pi.init
lambda.init <- init.values$lambda.init

## > pi.init
##  Cluster 1  Cluster 2  Cluster 3  Cluster 4 
## 0.09659369 0.33289211 0.38439010 0.18612411 
## >
## > lambda.init
##   Cluster 1 Cluster 2 Cluster 3 Cluster 4
## 1 0.2543011 0.7618805 1.2588717 1.2597545
## 2 0.9117585 1.5157468 0.9123846 0.3060517
## 3 2.8558488 0.5705673 0.5779289 1.7152452

}
\keyword{ models }

