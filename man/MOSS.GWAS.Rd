\name{MOSS.GWAS}
\Rdversion{1.0}
\title{A function implementing the MOSS algorithm for the analysis of GWAS data.}
\description{The MOSS algorithm is a Bayesian variable selection procedure that can be used for the analsysis GWAS data. It identifies combinations of the best predictive SNPs associated with the response. It also performs a hierarchical log-linear model search to identify the most relevant associations among the resulting subsets of SNPs. The function has an option to use model averaging to construct a classifier for predicting the response and to assess its capability via k-fold cross validation. The prior used is the generalized hyper Dirichlet.}
\usage{MOSS.GWAS(alpha = 1, c = 0.1, cPrime = 0.0001, q = 0.1,
          replicates = 5, maxVars = 3, data, dimens, k = 2)}
\alias{MOSS.GWAS}
\arguments{
  \item{alpha}{A hyperparameter of the prior representing the total of a fictive contingency table with counts equal to alpha divided by the number of cells. Alpha must be a positive real number.}
  \item{c, cPrime, q}{Tuning parameters for the MOSS algorithm. All 3 must be real numbers between 0 and 1 and cPrime must be smaller than c.} 
  \item{replicates}{The number of instances the MOSS algorithm will be run.}
  \item{maxVars}{The maximum number of variables allowed in a model (including the response). Must be an integer between 3 and 6.}
  \item{data}{A data frame containing the genotype information for a given set of SNPs. The data frame should be organized such that each row refers to a subject and each column to a SNP. The last column is interpreted as the case / control status of each subject and must be binary. Otherwise, the SNP data need not be binary.} 
  \item{dimens}{The number of possible values for each column of data. Each possible value does not need to occur in data. Since the last column of data must be binary, the last entry of dimens must be 2. All other entries of dimens must be greater than or equal to 2.} 
  \item{k}{The fold of the cross validation. If k is NULL then no cross validation is performed.}  
}
\value{A list with 4 data frame elements:
\item{topRegressions}{The top regressions identified together with their log marginal likelihood and normalized posterior probability.}
\item{postIncProbs}{The posterior inclusion probabilities of each variable that appears in one of the top regressions.}
\item{interactionModels}{The best (in terms of marginal likelihood) hierarchical log-linear model containing the variables in each of the top regressions.}
\item{crossValidation}{A table with the average results of the cross validation. This table is typically called a confusion matrix.}
}
\references{[1] Massam, H., Liu, J. and Dobra, A. (2009). A conjugate prior for discrete hierarchical log-linear models. Annals of Statistics, 37, 3431-3467.

  [2] Dobra, A., Briollais, L., Jarjanazi, H., Ozcelik, H. and Massam, H. (2010). Applications of the mode oriented stochastic search (MOSS) algorithm for discrete multi-way data to genomewide studies. Bayesian Modeling in Bioinformatics, Taylor & Francis (D. Dey, S. Ghosh and B. Mallick, eds.), 63-93. 

  [3] Dobra, A. and Massam, H. (2010). The mode oriented stochastic search (MOSS) algorithm for log-linear models with conjugate priors. Statistical Methodology, 7, 240-253. 
}
\author{Matthew Friedlander and Laurent Briollais} 
\examples{
data(simuCC)
data <- simuCC[,c(1002,2971,rep(5978:6001))]
# The SNPs in columns 1002 and 2971 of simuCC called rs4491689 and rs6869003 
# cause the disease
set.seed(7)
MOSS.GWAS (alpha = 1, c = 0.1, cPrime = 0.0001, q = 0.1, replicates = 1,
           maxVars = 3, data, dimens = c(rep(3,25),2), k = NULL)
}
\keyword{htest}
\keyword{models}
