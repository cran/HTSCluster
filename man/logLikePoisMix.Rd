\name{logLikePoisMix}
\alias{logLikePoisMix}
\alias{logLikePoisMixDiff}
\alias{mylogLikePoisMixObs}
\title{
Log likelihood calculation for a Poisson mixture model
}
\description{
Functions to calculate the log likelihood for a Poisson mixture model, the difference in log likelihoods for two different sets of parameters of a Poisson mixture model or the log-likelihood for each observation.
}
\usage{
logLikePoisMix(y, mean, pi)
logLikePoisMixDiff(y, mean.new, pi.new, mean.old, pi.old)
mylogLikePoisMixObs(y, conds, s, lambda, pi)
}
\arguments{
  \item{y}{
(\emph{n} x \emph{q}) matrix of observed counts for \emph{n} observations and \emph{q} variables
}
  \item{mean}{
List of length \emph{g} containing the (\emph{n} x \emph{q}) matrices of conditional mean expression for all observations, as calculated by the \code{\link{PoisMixMean}} function, where \emph{g} represents the number of clusters
}
  \item{mean.new}{
List of length \emph{g} containing the (\emph{n} x \emph{q}) matrices of conditional mean expression for all observations for one set of parameters, as calculated by the \code{\link{PoisMixMean}} function, where \emph{g} represents the number of clusters
}
  \item{mean.old}{
List of length \emph{g} containing the (\emph{n} x \emph{q}) matrices of conditional mean expression for all observations for another set of parameters, as calculated by the \code{\link{PoisMixMean}} function, where \emph{g} represents the number of clusters
}
  \item{pi.new}{
Vector of length \emph{g} containing one estimate for \eqn{\hat{\boldsymbol{\pi}}}{\hat{\pi}}
}
  \item{pi.old}{
Vector of length \emph{g} containing another estimate for \eqn{\hat{\boldsymbol{\pi}}}{\hat{\pi}}
}
  \item{pi}{
Vector of length \emph{g} containing estimate for \eqn{\hat{\boldsymbol{\pi}}}{\hat{\pi}}
}
  \item{conds}{
Vector of length \emph{q} defining the condition (treatment group) for each variable (column) in \code{y}
}
  \item{s}{
Estimate of normalized per-variable library size
}
  \item{lambda}{
(\emph{d} x \code{g}) matrix containing the current estimate of lambda, where \emph{d} is the number of conditions (treatment groups) and \code{g} is the number of clusters
}
}
\details{
The \code{logLikePoisMixDiff} function is used to calculate the difference in log likelihood for two different sets of parameters in a Poisson mixture model; it is used to determine convergence in the EM algorithm run by the \code{\link{PoisMixClus}} function.
The \code{logLikePoisMix} function (taken largely from the \code{mylogLikePoisMix} function from the \code{poisson.glm.mix} R package) calculates the log likelihood for a given set of parameters in a Poisson mixture model and is used in the \code{\link{PoisMixClus}} function for the calculation of the BIC and ICL. 
The \code{mylogLikePoisMixObs} function calculates the log likelihood per observation for a given set of parameters in a Poisson mixture model. 
}

\value{
\item{ll }{(Depending on the context), the log likelihood, difference in log likelihoods for two different sets of parameters, or per-observation log-likelihood}
}
\note{
In the \code{logLikePoisMixDiff} function, we make use of the alternative mass function for a Poisson density proposed by Loader (2000) to avoid computational difficulties. The \code{logLikePoisMixDiff} function returns a default value of 100 if one or both of the log likelihoods associated with the two parameter sets takes on a value of \eqn{-\infty}{-\infty}.
}
\references{

Loader, C. (2000) Fast and accurate computation of binomial probabilities. Available at \url{https://lists.gnu.org/archive/html/octave-maintainers/2011-09/pdfK0uKOST642.pdf}.

Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux G. (2015). Co-expression analysis of high-throughput transcriptome sequencing data with Poisson mixture models. Bioinformatics, 31(9):1420-1427.

Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011) Clustering high-throughput sequencing data with Poisson mixture models. Inria Research Report 7786. Available at \url{https://hal.inria.fr/inria-00638082}.
}
\author{
Andrea Rau
}

\seealso{
\code{\link{PoisMixClus}} for Poisson mixture model estimation and model selection; 
\code{\link{PoisMixMean}} to calculate the per-cluster conditional mean of each observation
}
\examples{

set.seed(12345)

## Simulate data as shown in Rau et al. (2011)
## Library size setting "A", low cluster separation
## n = 200 observations

simulate <- PoisMixSim(n = 200, libsize = "A", separation = "low")
y <- simulate$y
conds <- simulate$conditions
w <- rowSums(y)               ## Estimate of w
r <- table(conds)             ## Number of replicates per condition
d <- length(unique(conds))    ## Number of conditions
s <- colSums(y) / sum(y)      ## TC estimate of lib size
s.dot <- rep(NA, d)           ## Summing lib size within conditions
for(j in 1:d) s.dot[j] <- sum(s[which(conds == unique(conds)[j])]);

## Initial guess for pi and lambda
g.true <- 4
pi.guess <- simulate$pi
## Recalibrate so that (s.dot * lambda.guess) = 1
lambda.sim <- simulate$lambda
lambda.guess <- matrix(NA, nrow = d, ncol = g.true)
for(k in 1:g.true) {
    tmp <- lambda.sim[,k]/sum(lambda.sim[,k])
    lambda.guess[,k] <- tmp/s.dot
}

## Run the PMM-II model for g = 4
## with EM algorithm and "TC" library size parameter
run <- PoisMixClus(y, g = 4, norm = "TC", conds = conds) 
pi.est <- run$pi
lambda.est <- run$lambda

## Mean values for each of the parameter sets
mean.guess <- PoisMixMean(y, 4, conds, s, lambda.guess)
mean.est <- PoisMixMean(y, 4, conds, s, lambda.est)

## Difference in log likelihoods       
LL.diff <- logLikePoisMixDiff(y, mean.guess, pi.guess, mean.est, pi.est)
LL.diff             ## -12841.11

}
\keyword{ methods }

