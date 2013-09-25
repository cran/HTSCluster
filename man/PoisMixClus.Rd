\name{PoisMixClus}
\alias{PoisMixClus}
\alias{PoisMixClusWrapper}
\title{
Poisson mixture model estimation and model selection
}
\description{
These functions implement the EM and CEM algorithms for parameter estimation in a Poisson mixture model for clustering high throughput sequencing observations (e.g., genes) for a single number of clusters (\code{PoisMixClus}) or a sequence of cluster numbers (\code{PoisMixClusWrapper}). Parameters are initialized using a Small-EM strategy as described in Rau et al. (2011) or the splitting small-EM strategy described in Papastamoulis et al. (2012), and model selection is performed using the ICL criteria. Note that these functions implement the PMM-I and PMM-II models described in Rau et al. (2011). 
}
\usage{
PoisMixClus(y, g, conds, lib.size = TRUE, lib.type = "TMM", 
    init.type = "small-em", init.runs = 20, init.iter = 10,
    alg.type = "EM", cutoff = 10e-6, iter = 1000, fixed.lambda = NA,
    equal.proportions = FALSE, prev.labels = NA, 
    prev.probaPost = NA, verbose = FALSE, interpretation = "sum")

PoisMixClusWrapper(y, gmin = 1, gmax, conds, lib.size = TRUE,
    lib.type = "TMM", gmin.init.type = "small-em",
    init.runs = 20, init.iter = 10, split.init = TRUE, alg.type = "EM", 
    cutoff = 10e-6, iter = 1000, fixed.lambda = NA, 
    equal.proportions = FALSE, verbose = FALSE, interpretation = "sum")
}

\arguments{
  \item{y}{
(\emph{n} x \emph{q}) matrix of observed counts for \emph{n} observations and \emph{q} variables
}
  \item{g}{
Number of clusters (a single value). If \code{fixed.lambda} contains a list of lambda values to be fixed, 
\code{g} corresponds to the number of clusters in addition to those fixed.
}
  \item{gmin}{
The minimum number of clusters in a sequence to be tested, where \code{gmin} corresponds to the 
minimum number of clusters in addition to those fixed.
}
  \item{gmax}{
The maximum number of clusters in a sequence to be tested, where \code{gmax} corresponds to the 
maximum number of clusters in addition to those fixed.
}
  \item{conds}{
Vector of length \emph{q} defining the condition (treatment group) for each variable (column) in \code{y}
}
  \item{lib.size}{
If \code{FALSE}, the library size parameter is not included in the model (i.e., the PMM-I model). If \code{TRUE}, the library size parameter is included in the Poisson mixture model (i.e., the PMM-II model)
}
  \item{lib.type}{
If \code{lib.size = TRUE}, the estimator to be used for the library size parameter: (\dQuote{\code{TC}} for total count, \dQuote{\code{UQ}} for upper quantile, \dQuote{\code{Med}} for median, 
\dQuote{\code{DESeq}} for the normalization method in the DESeq package, and \dQuote{\code{TMM}} for the TMM normalization method.
}
  \item{init.type}{
Type of initialization strategy to be used (\dQuote{\code{small-em}} for the Small-EM strategy described in Rau et al. (2011), and \dQuote{\code{kmeans}} for a simple \emph{K}-means initialization)
}
  \item{gmin.init.type}{
Type of initialization strategy to be used for the minimum number of clusters in a sequence (\code{gmin}): 
(\dQuote{\code{small-em}} for the Small-EM strategy described in Rau et al. (2011), and \dQuote{\code{kmeans}} for a simple \emph{K}-means initialization)
}
  \item{init.runs}{
Number of iterations to be used for the Small-EM strategy described in Rau et al. (2011), with a default value of 20
}
  \item{init.iter}{
Number of iterations to be used within each run for the Small-EM strategry, with a default value of 10
}
  \item{split.init}{
If \code{TRUE}, the splitting initialization strategy of Papastamoulis et al. (2012) will be used for cluster sizes
(\code{gmin}+1, ..., \code{gmax}). If \code{FALSE}, the initialization strategy specified in \code{gmin.init.type}
is used for all cluster sizes in the sequence.
}
  \item{alg.type}{
Algorithm to be used for parameter estimation (\dQuote{\code{EM}} or \dQuote{\code{CEM}})
}
  \item{cutoff}{
Cutoff to declare algorithm convergence (in terms of differences in log likelihoods from one iteration to the next)
}
  \item{iter}{
Maximum number of iterations to be run for the chosen algorithm
}
  \item{fixed.lambda}{
If one (or more) clusters with fixed values of lambda is desired, a list containing vectors of length \emph{d} (the number of conditions). Note that the values of lambda chosen must satisfy the constraint noted in the technical report.
}
  \item{equal.proportions}{
If \code{TRUE}, the cluster proportions are set to be equal for all clusters. Default is \code{FALSE} (unequal cluster proportions).
}
  \item{prev.labels}{
A vector of length \emph{n} of cluster labels obtained from the previous run (g-1 clusters) to be used with the splitting
small-EM strategy described in described in Papastamoulis et al. (2012). For other initialization strategies, this 
parameter takes the value NA
}
  \item{prev.probaPost}{
An \emph{n} x (\emph{g}-1) matrix of the conditional probabilities of each observation belonging to each of the 
\emph{g}-1 clusters from the previous run, to be used with the splitting small-EM strategy of described in Papastamoulis 
et al. (2012). For other initialization strategies, this parameter takes the value NA
}
  \item{verbose}{
If \code{TRUE}, include verbose output
}
  \item{interpretation}{
If \code{"sum"}, cluster behavior is interpreted with respect to overall gene expression level (sums per gene), 
otherwise for \code{"mean"}, cluster behavior is interpreted with respect to mean gene expression (means per gene).
}
}
\details{
Output of \code{PoisMixClus} is an S3 object of class \code{HTSCluster}, and output of \code{PoisMixClusWrapper} is an S3 object
of class \code{HTSClusterWrapper}.

In a Poisson mixture model, the data \eqn{\mathbf{y}}{y} are assumed to come from \emph{g} distinct subpopulations (clusters), each of which is modeled separately; the overall population is thus a mixture of these subpopulations. In the case of a Poisson mixture model with \emph{g} components, the model may be written as

\deqn{f(\mathbf{y};g,\ensuremath\boldsymbol{\Psi}_g) = \prod_{i=1}^n \sum_{k=1}^g \pi_k \prod_{j=1}^{d}\prod_{l=1}^{r_j} P(y_{ijl} ; \ensuremath\boldsymbol{\theta}_k)}{f(y;g,\psi_g) = \prod_{i=1}^n \sum_{k=1}^g \pi_k \prod_{j=1}^{d}\prod_{l=1}^{r_j} P(y_{ijl} ; \theta_k)}

for \eqn{i = 1, \ldots, n} observations in \eqn{l = 1, \ldots, r_j} replicates of \eqn{j = 1, \ldots, d} conditions (treatment groups), where \eqn{P(\cdot)} is the standard Poisson density, \eqn{\ensuremath\boldsymbol{\Psi}_g = (\pi_1,\ldots,\pi_{g-1}, \ensuremath\boldsymbol{\theta}^\prime)}{\psi_g = (\pi_1,\ldots,\pi_{g-1}, \theta^\prime)}, \eqn{\ensuremath\boldsymbol{\theta}^\prime}{\theta^\prime} contains all of the parameters in \eqn{\ensuremath\boldsymbol{\theta}_1,\ldots,\ensuremath\boldsymbol{\theta}_g}{\theta_1,\ldots,\theta_g} assumed to be distinct, and \eqn{\ensuremath\boldsymbol{\pi} = (\pi_1,\ldots,\pi_g)^\prime}{\pi = (\pi_1,\ldots,\pi_g)^\prime} are the mixing proportions such that \eqn{\pi_k} is in (0,1) for all \emph{k} and \eqn{\sum_k \pi_k = 1}.

We consider two possible parameterizations for the mean \eqn{\ensuremath\boldsymbol{\theta}_k = (\mu_{ijlk})}{\theta = (mu_{ijlk})}. In the first, called the PMM-I, we consider
\deqn{\mu_{ijlk} = w_i \lambda_{jk}}
where \eqn{w_i} corresponds to the expression level of observation \emph{i} and \eqn{\ensuremath\boldsymbol{\lambda}_k = (\lambda_{1k},\ldots,\lambda_{dk})}{\lambda_k = (\lambda_{1k},\ldots,\lambda_{dk})} corresponds to the clustering parameters that define the profiles of the genes in cluster \emph{k} across all variables. In the second parameterization, called the PMM-II, we consider
\deqn{\mu_{ijlk} = w_i s_{jl} \lambda_{jk}}
where \eqn{w_i} and \eqn{\ensuremath\boldsymbol{\lambda}_k}{\lambda_k} are as before and \eqn{s_{jl}} is the normalized library size (a fixed constant) for replicate \emph{l} of condition \emph{j}. See Rau et al. (2011) for more details on the PMM-I and PMM-II.

There are two approaches to estimating the parameters of a finite mixture model and obtaining a clustering of the data: the estimation approach (via the EM algorithm) and the clustering approach (via the CEM algorithm). Parameter initialization is done using a Small-EM strategy as described in Rau et al. (2011) via the \code{\link{emInit}} function. Model selection may be performed using the BIC or ICL criteria.
}
\value{
\item{lambda }{(\emph{d} x \emph{g}) matrix containing the estimate of \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}}}
\item{pi }{Vector of length \emph{g} containing the estimate of \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}}}
\item{labels }{Vector of length \emph{n} containing the cluster assignments of the \emph{n} observations}
\item{probaPost }{Matrix containing the conditional probabilities of belonging to each cluster for all observations}
\item{log.like }{Value of log likelihood}
\item{BIC }{Value of BIC criterion}
\item{ICL }{Value of ICL criterion}
\item{alg.type }{Estimation algorithm used; matches the argument \code{alg.type} above)}
\item{lib.size }{TRUE if library size included in the model (matches the argument \code{alg.type} above)}
\item{lib.type }{Type of library size normalization used (if \code{lib.size = TRUE}; matches the argument \code{lib.type} above)}
\item{s }{Library size normalization factors used (if \code{lib.size = TRUE})}
\item{conds }{Conditions specified by user}
\item{loglike.all }{Log likelihoods calculated for each of the fitted models for cluster sizes \code{gmin}, ..., \code{gmax}}
\item{ICL.all }{ICL values calculated for each of the fitted models for cluster sizes \code{gmin}, ..., \code{gmax}}
\item{select.results }{Object of class \code{HTSCluster} giving the results from the model chosen via the ICL criterion}
\item{all.results }{List of objects of class \code{HTSCluster} giving the results for all models for cluster sizes \code{gmin}, ..., \code{gmax}}
}

\references{
Anders, S. and Huber, W. (2010) Differential expression analysis for sequence count data. \emph{Genome Biology}, \bold{11}(R106), 1-28.

Papastamoulis, P., Martin-Magniette, M.-L., and Maugis-Rabusseau, C. (2012). Efficient estimation of high dimensional mixtures of Poisson generalized linear models via the EM algorithm (submitted).

Rau, A., Celeux, G., Martin-Magniette, M.-L., Maugis-Rabusseau, C (2011). Clustering high-throughput sequencing data with Poisson mixture models. Inria Research Report 7786. Available at \url{http://hal.inria.fr/inria-00638082}.
}
\author{
Andrea Rau <\url{andrea.rau@jouy.inra.fr}>
}

\seealso{
\code{\link{probaPost}} for the calculation of the conditional probability of belonging to a cluster;
\code{\link{PoisMixMean}} for the calculation of the per-cluster conditional mean of each observation;
\code{\link{logLikePoisMix}} for the calculation of the log likelihood of a Poisson mixture model;
\code{\link{emInit}} and \code{\link{kmeanInit}} for the Small-EM parameter initialization strategy
}
\examples{

set.seed(12345)

## Simulate data as shown in Rau et al. (2011)
## Library size setting "A", high cluster separation
## n = 200 observations

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
\keyword{ models }
\keyword{ cluster }

