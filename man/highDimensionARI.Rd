\name{highDimensionARI}
\alias{highDimensionARI}
\title{
Calculate ARI for high-dimensional data via data splits
}
\description{
This function is used to calculate Adjusted Rand Index (ARI) values for high-dimensional data.
}
\usage{
highDimensionARI(x, y, splits = 2, verbose = FALSE) 
}
\arguments{
  \item{x}{
Vector of classification labels
}
  \item{y}{
Vector of classification labels
}
  \item{splits}{
Number of subsets data should be split into
}
  \item{verbose}{
\code{TRUE} if verbose output is desired
}
}

\value{
Value of Adjusted Rand Index for samples \code{x} and \code{y}
}
\references{
Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011). Clustering high-throughput sequencing data with Poisson mixture models. Inria Research Report 7786. Available at \url{http://hal.inria.fr/inria-00638082}.
}
\author{
Andrea Rau <\url{andrea.rau@jouy.inra.fr}>
}
\keyword{ methods }

