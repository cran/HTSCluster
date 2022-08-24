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
Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux G. (2015). Co-expression analysis of high-throughput transcriptome sequencing data with Poisson mixture models. Bioinformatics, 31(9):1420-1427.

Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011). Clustering high-throughput sequencing data with Poisson mixture models. Inria Research Report 7786. Available at \url{https://hal.inria.fr/inria-00638082}.
}
\author{
Andrea Rau
}
\keyword{ methods }

