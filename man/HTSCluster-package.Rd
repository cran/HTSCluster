\name{HTSCluster-package}
\alias{HTSCluster-package}
\alias{HTSCluster}
\docType{package}
\title{
Clustering high throughput sequencing (HTS) data
}
\description{
A Poisson mixture model is implemented to cluster genes from high-throughput 
transcriptome sequencing (RNA-seq) data. Parameter estimation is performed using 
either the EM or CEM algorithm, and the slope heuristics are used for model 
selection (i.e., to choose the number of clusters).
}
\author{
Andrea Rau, Gilles Celeux, Marie-Laure Martin-Magniette, Cathy Maugis-Rabusseau

Maintainer: Andrea Rau
}
\references{
Rau, A., Maugis-Rabusseau, C., Martin-Magniette, M.-L., Celeux G. (2015). Co-expression analysis of high-throughput transcriptome sequencing data with Poisson mixture models. Bioinformatics, 31(9):1420-1427.

Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C. (2011) Clustering high-throughput sequencing data with Poisson mixture models. Inria Research Report 7786. Available at \url{https://inria.hal.science/inria-00638082}.
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

run <- PoisMixClus(y, g=3, conds=conds, norm="TC") 

## Estimates of pi and lambda for the selected model

pi.est <- run$pi
lambda.est <- run$lambda


## Not run: PMM for 4 total clusters, with one fixed class
## "TC" library size estimate, EM algorithm
##
## run <- PoisMixClus(y, g = 3, norm = "TC", conds = conds, 
##    fixed.lambda = list(c(1,1,1))) 
##
##
## Not run: PMM model for 4 clusters, with equal proportions
## "TC" library size estimate, EM algorithm
##
## run <- PoisMixClus(y, g = 4, norm = "TC", conds = conds, 
##     equal.proportions = TRUE) 
##
##
## Not run: PMM model for g = 1, ..., 10 clusters, Split Small-EM init
##
## run1.10 <- PoisMixClusWrapper(y, gmin = 1, gmax = 10, conds = conds, 
##	norm = "TC")
##
##
## Not run: PMM model for g = 1, ..., 10 clusters, Small-EM init
##
## run1.10bis <-  <- PoisMixClusWrapper(y, gmin = 1, gmax = 10, conds = conds, 
##	norm = "TC", split.init = FALSE)
##
##
## Not run: previous model equivalent to the following
##
## for(K in 1:10) {
##	run <- PoisMixClus(y, g = K, conds = conds, norm = "TC")
## }


}
