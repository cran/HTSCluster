\name{probaPost}
\alias{probaPost}
\title{
Calculate the conditional probability of belonging to each cluster in a Poisson mixture model
}
\description{
This function computes the conditional probabilities \eqn{t_{ik}} that an observation \emph{i} arises from the \eqn{k^{\mathrm{th}}}{kth} component for the current value of the mixture parameters.
}
\usage{
probaPost(y, g, conds, pi, s, lambda)
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
  \item{pi}{
Vector of length \code{g} containing the current estimate of \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\pi}
}
  \item{s}{
Vector of length \emph{q} containing the estimates for the normalized library size parameters for each of the \emph{q} variables in \code{y}
}
  \item{lambda}{
(\emph{d} x \code{g}) matrix containing the current estimate \eqn{\ensuremath\boldsymbol{\lambda}}{\lambda}, where \emph{d} is the number of conditions (treatment groups)
}
}

\value{
\item{t }{(\emph{n} x \code{g}) matrix made up of the conditional probability of each observation belonging to each of the \code{g} clusters}

}
\references{
Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux G. (2015). Co-expression analysis of high-throughput transcriptome sequencing data with Poisson mixture models. Bioinformatics, 31(9):1420-1427.

Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011). Clustering high-throughput sequencing data with Poisson mixture models. Inria Research Report 7786. Available at \url{http://hal.inria.fr/inria-00638082}.
}
\author{
Andrea Rau <\url{andrea.rau@jouy.inra.fr}>
}
\note{
If all values of \eqn{t_{ik}} are 0 (or nearly zero), the observation is assigned with probability one to belong to the cluster with the closest mean (in terms of the Euclidean distance from the observation). To avoid calculation difficulties, extreme values of \eqn{t_{ik}} are smoothed, such that those smaller than 1e-10 or larger than 1-1e-10 are set equal to 1e-10 and 1-1e-10, respectively.
}


\seealso{
\code{\link{PoisMixClus}} for Poisson mixture model estimation and model selection;
\code{\link{PoisMixMean}} to calculate the conditional per-cluster mean of each observation
}
\examples{

set.seed(12345)

## Simulate data as shown in Rau et al. (2011)
## Library size setting "A", high cluster separation
## n = 200 observations

simulate <- PoisMixSim(n = 200, libsize = "A", separation = "high")
y <- simulate$y
conds <- simulate$conditions
s <- colSums(y) / sum(y)     ## TC estimate of lib size

## Run the PMM-II model for g = 3
## "TC" library size estimate, EM algorithm

run <- PoisMixClus(y, g = 3, norm = "TC",
 	conds = conds) 
pi.est <- run$pi
lambda.est <- run$lambda

## Calculate the conditional probability of belonging to each cluster
proba <- probaPost(y, g = 3, conds = conds, pi = pi.est, s = s, 
	lambda = lambda.est)

## head(round(proba,2))

}
\keyword{ methods }

