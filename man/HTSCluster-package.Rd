\name{HTSCluster-package}
\alias{HTSCluster-package}
\alias{HTSCluster}
\docType{package}
\title{
Clustering high throughput sequencing (HTS) data
}
\description{
This package implements estimation of a Poisson mixture model to 
cluster observations (e.g., genes) in high throughput sequencing data. 
Parameter estimation is performed using either the EM or CEM algorithm, and the 
slope heuristics are used for model selection (i.e., to choose the number of 
clusters).
}
\details{
\tabular{ll}{
Package: \tab HTSCluster\cr
Type: \tab Package\cr
Version: \tab 2.0.4\cr
Date: \tab 2014-08-26\cr
License: \tab GPL (>=3)\cr
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

## Run the PMM model for g = 3
## "TC" library size estimate, EM algorithm

run <- PoisMixClus(y, g = 3, conds = conds, lib.type = "TC") 

## Estimates of pi and lambda for the selected model

pi.est <- run$pi
lambda.est <- run$lambda


## Not run: PMM for 4 total clusters, with one fixed class
## "TC" library size estimate, EM algorithm
##
## run <- PoisMixClus(y, g = 3, lib.type = "TC", conds = conds, 
##    fixed.lambda = list(c(1,1,1))) 
##
##
## Not run: PMM model for 4 clusters, with equal proportions
## "TC" library size estimate, EM algorithm
##
## run <- PoisMixClus(y, g = 4, lib.type = "TC", conds = conds, 
##     equal.proportions = TRUE) 
##
##
## Not run: PMM model for g = 1, ..., 10 clusters, Split Small-EM init
##
## run1.10 <- PoisMixClusWrapper(y, gmin = 1, gmax = 10, conds = conds, 
##	lib.type = "TC")
##
##
## Not run: PMM model for g = 1, ..., 10 clusters, Small-EM init
##
## run1.10bis <-  <- PoisMixClusWrapper(y, gmin = 1, gmax = 10, conds = conds, 
##	lib.type = "TC", split.init = FALSE)
##
##
## Not run: previous model equivalent to the following
##
## for(K in 1:10) {
##	run <- PoisMixClus(y, g = K, conds = conds, lib.type = "TC")
## }


}
