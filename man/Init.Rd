\name{Init}
\alias{emInit}
\alias{kmeanInit}
\alias{splitEMInit}
\alias{probaPostInit}
\title{
Parameter initialization for a Poisson mixture model.
}
\description{
These functions implement a variety of initialization methods for the parameters of a Poisson mixture model: the Small EM initialization strategy (\code{emInit}) described in Rau et al. (2011), a K-means initialization strategy (\code{kmeanInit}) that is itself used to initialize the small EM strategy, the splitting small-EM initialization strategy (\code{splitEMInit}) based on that described in Papastamoulis et al. (2014), and a function to initialize a small-EM strategy using the posterior probabilities (\code{probaPostInit}) obtained from a previous run with one fewer cluster following the splitting strategy.
}
\usage{

emInit(y, g, conds, norm, alg.type = "EM", 
    init.runs, init.iter, fixed.lambda, equal.proportions, verbose)

kmeanInit(y, g, conds, norm, fixed.lambda,
    equal.proportions)

splitEMInit(y, g, conds, norm, alg.type, fixed.lambda,
    equal.proportions, prev.labels, prev.probaPost, init.runs, 
    init.iter, verbose)

probaPostInit(y, g, conds, norm, alg.type = "EM", 
    fixed.lambda, equal.proportions, probaPost.init, init.iter,
    verbose)

}
\arguments{
  \item{y}{
(\emph{n} x \emph{q}) matrix of observed counts for \emph{n} observations and \emph{q} variables
}
  \item{g}{
Number of clusters. If \code{fixed.lambda} contains a list of lambda values to be fixed, \emph{g} corresponds to the number of clusters in addition to those fixed.
}
  \item{conds}{
Vector of length \emph{q} defining the condition (treatment group) for each variable (column) in \code{y}
}
  \item{norm}{
The type of estimator to be used to normalize for differences in library size: (\dQuote{\code{TC}} for total count, \dQuote{\code{UQ}} for upper quantile, \dQuote{\code{Med}} for median, 
\dQuote{\code{DESeq}} for the normalization method in the DESeq package, and \dQuote{\code{TMM}} for the TMM normalization method (Robinson and Oshlack, 2010). Can also
be a vector (of length \emph{q}) containing pre-estimated library size estimates for each sample.
}
  \item{alg.type}{
Algorithm to be used for parameter estimation (\dQuote{\code{EM}} or \dQuote{\code{CEM}} for the EM or CEM algorithms, respectively)
}
  \item{init.runs}{
In the case of the Small-EM algorithm, the number of independent runs to be performed.
In the case of the splitting Small-EM algorithm, the number of cluster splits to be performed in the splitting small-EM initialization. 
}
  \item{init.iter}{
The number of iterations to run within each Small-EM algorithm
}
  \item{fixed.lambda}{
If one (or more) clusters with fixed values of lambda is desires, a list containing vectors of length \emph{d} (the number of conditions). Note that the values of lambda chosen must satisfy the constraint noted in the technical report.
}
  \item{equal.proportions}{
If \code{TRUE}, the cluster proportions are set to be equal for all clusters. Default is \code{FALSE} (unequal cluster proportions)
}
  \item{prev.labels}{
A vector of length \emph{n} of cluster labels obtained from the previous run (g-1 clusters)
}
  \item{prev.probaPost}{
An \emph{n} x (\emph{g}-1) matrix of the conditional probabilities of each observation belonging to each of the 
\emph{g}-1 clusters from the previous run
}
  \item{probaPost.init}{
An \emph{n} x (\emph{g}) matrix of the conditional probabilities of each observation belonging to each of the 
\emph{g} clusters following the splitting strategy in the \code{splitEMInit} function
}
  \item{verbose}{
If \code{TRUE}, include verbose output
}
}
\details{
In practice, the user will not directly call the initialization functions described here; they are indirectly called
for a single number of clusters through the \code{PoisMixClus} function (via \code{init.type}) or via the 
\code{PoisMixClusWrapper} function for a sequence of cluster numbers (via \code{gmin.init.type} and \code{split.init}). 

To initialize parameter values for the EM and CEM algorithms, for the Small-EM strategy (Biernacki et al., 2003) we use the \code{emInit} function as follows. For a given number of independent runs (given by \code{init.runs}), the following procedure is used to obtain parameter values: first, a K-means algorithm (MacQueen, 1967) is run to partition the data into \code{g} clusters (\eqn{\hat{\ensuremath\boldsymbol{z}}^{(0)}}{\hat{z}^(0)}). Second, initial parameter values \eqn{\ensuremath\boldsymbol{\pi}^{(0)}}{\pi^(0)} and \eqn{\ensuremath\boldsymbol{\lambda}^{(0)}}{\lambda^(0)} are calculated (see Rau et al. (2011) for details). Third, a given number of iterations of an EM algorithm are run (defined by \code{init.iter}), using \eqn{\ensuremath\boldsymbol{\pi}^{(0)}}{\pi^(0)} and \eqn{\ensuremath\boldsymbol{\lambda}^{(0)}}{\lambda^(0)} as initial values. Finally, among the \code{init.runs} sets of parameter values, we use \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}} and \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}} corresponding to the highest log likelihood or completed log likelihood to initialize the subsequent full EM or CEM algorithms, respectively.

For the splitting small EM initialization strategy, we implement an approach similar to that described in Papastamoulis et al. (2014),
where the cluster from the previous run (with \emph{g}-1 clusters) with the largest entropy is chosen to be split into two new clusters,
followed by a small EM run as described above.
}
\value{
\item{pi.init }{Vector of length \code{g} containing the estimate for \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}} corresponding to the highest log likelihood (or completed log likelihood) from the chosen inialization strategy. }

\item{lambda.init }{(\emph{d} x \code{g}) matrix containing the estimate of \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}} corresponding to the highest log likelihood (or completed log likelihood) from the chosen initialization strategy, where \emph{d} is the number of conditions and \code{g} is the number of clusters. }

\item{lambda }{(\emph{d} x \code{g}) matrix containing the estimate of \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}} arising from the splitting initialization and small EM run for a single split, 
where \emph{d} is the number of conditions and \code{g} is the number of clusters. }

\item{pi }{Vector of length \code{g} containing the estimate for \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}} arising from the splitting initialization and small EM run for a single split, where \code{g} is the 
number of clusters.  }

\item{log.like }{Log likelihood arising from the splitting initialization and small EM run for a single split. }

}
\references{
Anders, S. and Huber, W. (2010) Differential expression analysis for sequence count data. \emph{Genome Biology}, \bold{11}(R106), 1-28.

Biernacki, C., Celeux, G., Govaert, G. (2003) Choosing starting values for the EM algorithm for getting the highest likelhiood in multivariate Gaussian mixture models. \emph{Computational Statistics and Data Analysis}, \bold{41}(1), 561-575.

MacQueen, J. B. (1967) Some methods for classification and analysis of multivariate observations. In \emph{Proceedings of the 5th Berkeley Symposium on Mathematical Statistics and Probability}, number 1, pages 281-297. Berkeley, University of California Press.

Papastamoulis, P., Martin-Magniette, M.-L., and Maugis-Rabusseau, C. (2014). On the estimation of mixtures of Poisson regression models with large number of components. \emph{Computational Statistics and Data Analysis}: 3rd special Issue on Advances in Mixture Models, DOI: 10.1016/j.csda.2014.07.005.

Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011). Clustering high-throughput sequencing data with Poisson mixture models. Inria Research Report 7786. Available at \url{http://hal.inria.fr/inria-00638082}.

Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux G. (2015). Co-expression analysis of high-throughput transcriptome sequencing data with Poisson mixture models. Bioinformatics, 31(9):1420-1427.

Robinson, M. D. and Oshlack, A. (2010) A scaling normalization method for differential expression analysis of RNA-seq data. \emph{Genome Biology}, \bold{11}(R25).
}
\author{
Andrea Rau <\url{andrea.rau@jouy.inra.fr}>
}

\seealso{
\code{\link{PoisMixClus}} for Poisson mixture model estimation for a given number of clusters,
\code{\link{PoisMixClusWrapper}} for Poisson mixture model estimation and model selection for a sequence of cluster numbers.
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
##
## init.values <- emInit(y, g = 4, conds, 
##    norm = "TC", alg.type = "EM", 
##    init.runs = 50, init.iter = 10, fixed.lambda = NA,
##    equal.proportions = FALSE, verbose = FALSE)
## pi.init <- init.values$pi.init
## lambda.init <- init.values$lambda.init


}
\keyword{ models }

