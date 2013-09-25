\name{PoisMixMean}
\alias{PoisMixMean}
\title{
Calculate the conditional per-cluster mean of each observation
}
\description{
This function is used to calculate the conditional per-cluster mean expression for all observations. This value corresponds to
\eqn{\ensuremath\boldsymbol{\mu} = (\mu_{ijlk}) = (\hat{w}_i \hat{\lambda}_{jk})}{\mu = (\mu_{ijlk}) = (\hat{w}_i \hat{\lambda}_{jk})}
for the PMM-I model and
\eqn{\ensuremath\boldsymbol{\mu} = (\mu_{ijlk}) = (\hat{w}_i s_{jl} \hat{\lambda}_{jk})}{\mu = (\mu_{ijlk}) = (\hat{w}_i s_{jl}\hat{\lambda}_{jk})}
for the PMM-II model.
}
\usage{
PoisMixMean(y, g, conds, s, lambda)
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
  \item{s}{
Estimate of normalized per-variable library size
}
  \item{lambda}{
(\emph{d} x \code{g}) matrix containing the current estimate of lambda, where \emph{d} is the number of conditions (treatment groups) and \code{g} is the number of clusters
}
}

\value{
A list of length \code{g} containing the (\emph{n} x \emph{q}) matrices of mean expression for all observations, conditioned on each of the \code{g} clusters
}
\references{
Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011). Clustering high-throughput sequencing data with Poisson mixture models. Inria Research Report 7786. Available at \url{http://hal.inria.fr/inria-00638082}.
}
\author{
Andrea Rau <\url{andrea.rau@jouy.inra.fr}>
}
\seealso{
\code{\link{PoisMixClus}} for Poisson mixture model estimation and model selection
}
\examples{

set.seed(12345)

## Simulate data as shown in Rau et al. (2011)
## Library size setting "A", high cluster separation
## n = 200 observations

simulate <- PoisMixSim(n = 200, libsize = "A", separation = "high")
y <- simulate$y
conds <- simulate$conditions
s <- colSums(y) / sum(y) 	## TC estimate of lib size

## Run the PMM-II model for g = 3
## "TC" library size estimate, EM algorithm

run <- PoisMixClus(y, g = 3, lib.size = TRUE, 
    lib.type = "TC", conds = conds) 
pi.est <- run$pi
lambda.est <- run$lambda

## Calculate the per-cluster mean for each observation
means <- PoisMixMean(y, g = 3, conds, s, lambda.est) 

}
\keyword{ methods }

