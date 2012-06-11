\name{PoisMixClus}
\alias{PoisMixClus}
\title{
Poisson mixture model estimation and model selection
}
\description{
This function implements the EM and CEM algorithms for parameter estimation in a Poisson mixture model for clustering high throughput sequencing observations (e.g., genes). Parameters are initialized using a Small-EM strategy as described in Rau et al. (2011), and model selection may be performed using the BIC and ICL criteria. Note that this function implements the PMM-I and PMM-II models described in Rau et al. (2011). 
}
\usage{
PoisMixClus(y, gmin, gmax, lib.size = TRUE, lib.type = "TC", 
    conds, init.type = "small-em", alg.type = "EM", 
    cutoff = 1e-05, iter = 1000, mean.filter = FALSE, verbose = FALSE)
}

\arguments{
  \item{y}{
(\emph{n} x \emph{q}) matrix of observed counts for \emph{n} observations and \emph{q} variables
}
  \item{gmin}{
Minimum number of clusters to be fit (must be less than or equal to \code{gmax})
}
  \item{gmax}{
Maximum number of clusters to be fit (must be greater than or equal to \code{gmax})
}
  \item{lib.size}{
If \code{FALSE}, the library size parameter is not included in the model (i.e., the PMM-I model). If \code{TRUE}, the library size parameter is included in the Poisson mixture model (i.e., the PMM-II model)
}
  \item{lib.type}{
If \code{lib.size = TRUE}, the estimator to be used for the library size parameter (\dQuote{\code{TC}} for total count, \dQuote{\code{Q}} for quantile, and \dQuote{\code{MedRatio}} for the median ratio of Anders and Huber (2010))
}
  \item{conds}{
Vector of length \emph{q} defining the condition (treatment group) for each variable (column) in \code{y}
}
  \item{init.type}{
Type of initialization strategy to be used (\dQuote{\code{small-em}} for the Small-EM strategy described in Rau et al. (2011), and \dQuote{\code{kmeans}} for a simple \emph{K}-means initialization)
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
  \item{mean.filter}{
Option threshold value for filtering genes with a mean count across all samples less than \code{mean.filter}.
}
  \item{verbose}{
If \code{TRUE}, include verbose output
}
}
\details{
Output is an S3 object of class \code{HTSCluster}.

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
\item{lambda }{List of length (\code{gmax}-\code{gmin}+1) containing the estimates \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}} for each of the models (\emph{g} = \code{gmin}, \ldots, \code{gmax}). For a given model \emph{g}, \code{lambda[[g]]} is a matrix of dimension (\emph{d} x \emph{g}), where \emph{d} is the number of conditions (treatment groups) and \emph{g} is the number of clusters}
\item{pi }{List of length (\code{gmax}-\code{gmin}+1) containing the estimates \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}} for each of the models (\emph{g} = \code{gmin}, \ldots, \code{gmax}). For a given model \emph{g}, \code{pi[[g]]} is a vector of length \emph{g}, where \emph{g} is the number of clusters}
\item{labels }{Matrix of dimension (\emph{n} x (\code{gmax}-\code{gmin}+1)) containing the label assignments for each of the \emph{n} observations in each of the models (\emph{g} = \code{gmin}, \ldots, \code{gmax})}
\item{probaPost }{List of length (\code{gmax}-\code{gmin}+1) containing the conditional probabilities of belonging to each cluster for all observations for each of the models (\emph{g} = \code{gmin}, \ldots, \code{gmax}). For a given model \emph{g}, \code{probaPost[[g]]} is a (\emph{n} x \emph{g}) matrix containing the conditional probailities of belonging to each of the \emph{g} clusters for all observations}
\item{BIC.all }{BIC values for each of the models (\emph{g} = \code{gmin}, \ldots, \code{gmax})}
\item{ICL.all }{ICL values for each of the models (\emph{g} = \code{gmin}, \ldots, \code{gmax})}
\item{alg.type }{Algorithm used for parameter estimation (matches the argument \code{alg.type} above)}
\item{BIC }{Maximum BIC value across all models considered (\emph{g} = \code{gmin}, \ldots, \code{gmax})}
\item{ICL }{Maximum ICL value across all models considered (\emph{g} = \code{gmin}, \ldots, \code{gmax})}
\item{g.BIC }{Number of clusters corresponding to the maximum BIC value}
\item{g.ICL }{Number of clusters corresponding to the maximum ICL value}
\item{labels.BIC }{Vector of length \emph{n} containing the cluster assignments of the \emph{n} observations in the model selected via the BIC (with number of clusters \code{g.BIC})}
\item{labels.ICL }{Vector of length \emph{n} containing the cluster assignments of the \emph{n} observations in the model selected via the ICL (with number of clusters \code{g.ICL})}
\item{lambda.BIC }{(\emph{d} x \code{g.BIC}) matrix containing the estimate of \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}} for the model selected via the BIC (with number of clusters \code{g.BIC})}
\item{pi.BIC }{Vector of length \code{g.BIC} containing the estimate of \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}} for the model selected via the BIC (with number of clusters \code{g.BIC})}
\item{lambda.ICL }{(\emph{d} x \code{g.ICL}) matrix containing the estimate of \eqn{\hat{\ensuremath\boldsymbol{\lambda}}}{\hat{\lambda}} for the model selected via the ICL (with number of clusters \code{g.ICL})}
\item{pi.ICL }{Vector of length \code{g.ICL} containing the estimate of \eqn{\hat{\ensuremath\boldsymbol{\pi}}}{\hat{\pi}} for the model selected via the ICL (with number of clusters \code{g.ICL})}
\item{probaPost.BIC }{(\emph{n} x \code{g.BIC}) matrix containing the conditional probabilities of belonging to each cluster for all observations for the model selected via the BIC (with number of clusters \code{g.BIC})}
\item{probaPost.ICL }{(\emph{n} x \code{g.ICL}) matrix containing the conditional probabilities of belonging to each cluster for all observations for the model selected via the ICL (with number of clusters \code{g.ICL})}
\item{lib.size }{TRUE if library size included in the model (matches the argument \code{alg.type} above)}
\item{lib.type }{Type of library size normalization used (if \code{lib.size = TRUE}; matches the argument \code{alg.type} above)}
\item{s }{Library size normalization factors used (if \code{lib.size = TRUE})}
\item{y }{Data (matches the argument \code{y} above, unless observations were removed using the optional \code{mean.filter})}
\item{ind.remove }{Indices for observations removed using the optional \code{mean.filter}}
\item{mean.filter }{Matched the argument \code{mean.filter} above}
}

\references{
Anders, S. and Huber, W. (2010) Differential expression analysis for sequence count data. \emph{Genome Biology}, \bold{11}(R106), 1-28.

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

## Run the PMM-II model for g = {3, 4, 5}
## "TC" library size estimate, EM algorithm
## Model selection via the ICL

run <- PoisMixClus(y, gmin = 3, gmax = 5, lib.size = TRUE, lib.type = "TC",
    conds = conds, init.type = "small-em") 

## Estimates of pi and lambda for the selected model
pi.est <- run$pi.ICL
lambda.est <- run$lambda.ICL

}
\keyword{ models }
\keyword{ cluster }

