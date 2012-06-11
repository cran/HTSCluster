\name{HTSCluster-package}
\alias{HTSCluster-package}
\alias{HTSCluster}
\docType{package}
\title{
Clustering high throughput sequencing (HTS) data
}
\description{
This package implements two parameterizations of a Poisson mixture model to 
cluster observations (e.g., genes) in high throughput sequencing data. 
Parameter estimation is performed using either the EM or CEM algorithm, and the 
BIC or ICL criteria are used for model selection (i.e., to choose the number of 
clusters).
}
\details{
\tabular{ll}{
Package: \tab HTSCluster\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2012-06-11\cr
License: \tab GPM (>=2)\cr
LazyLoad: \tab yes\cr
}
}
\author{
Andrea Rau, Gilles Celeux, Marie-Laure Martin-Magniette, Cathy Maugis-Rabusseau

Maintainer: Andrea Rau <\url{andrea.rau@jouy.inra.fr}>
}
\references{
Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C (2011). Clustering high-throughput sequencing data with Poisson mixture models. Inria Research Report 7786. Available at \url{http://hal.inria.fr/inria-00638082}.
}
\keyword{ models }
\keyword{ cluster }
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
## Model selection via the ICL

run <- PoisMixClus(y, gmin = 3, gmax = 5, lib.size = TRUE, 
    lib.type = "TC", conds = conds, init.type = "small-em") 

## Estimates of pi and lambda for the selected model
pi.est <- run$pi.ICL
lambda.est <- run$lambda.ICL

}
