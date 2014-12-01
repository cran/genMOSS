\name{MOSS_GWAS}
\Rdversion{1.2}
\title{Analyzing GWAS data using MOSS}
\description{The Mode Oriented Stochastic Search (MOSS) of Dobra at. al (2010) is a two stage Bayesian variable selection procedure for analyzing GWAS data. If we let Y be a response, X be a small group of SNPs, and r = Y|X be called the regression of Y on X, the first stage of the procedure is to find regressions such that the marginal likelihood, P(r) = P(Y|X), is high. In the second stage, a hierarchical log-linear model search is performed to identify the most relevant interactions among the variables in each of the top regressions. Using the top log-linear models, model averaging is then used to construct a classifier for predicting the response and its capability is assessed via k-fold cross validation. The prior used in Bayesian computations is the generalized hyper Dirichlet of Massam et. al (2009).}
\usage{MOSS_GWAS(alpha = 1, c = 0.1, cPrime = 0.0001, q = 0.1,
          replicates = 5, maxVars = 3, data, dimens, confVars = NULL, k = NULL)}
\alias{MOSS_GWAS}
\arguments{
  \item{alpha}{A hyperparameter of the prior representing the total of a fictive contingency table with counts equal to alpha divided by the number of cells. Alpha must be a positive real number.}
  \item{c, cPrime, q}{Tuning parameters for MOSS. All 3 must be real numbers between 0 and 1 and cPrime must be smaller than c.} 
  \item{replicates}{The number of instances the first stage of the MOSS procedure will be run. The top regressions are culled from the results of all the replicates.}
  \item{maxVars}{The maximum number of variables allowed in a regression (including the response). Must be an integer from 3 to 6.}
  \item{data}{A data frame containing the genotype information for a given set of SNPs. The data frame should be organized such that each row refers to a subject and each column to a SNP. The last column must be a binary response for each subject. The data frame must contain at least 8 columns. Rows containing any missing values (i.e. NAs) are omitted from the analysis.} 
  \item{dimens}{The number of possible values for each column of data. Each possible value does not need to occur in data. All entries of dimens must be greater than or equal to 2.} 
  \item{confVars}{The parameter confVars (for confounding variables) is a character vector specifying the names of SNPs which, other than the response, will be forced to be in every regression. If no confounding variables are desired, confVars can be set to NULL. A maximum of (maxVars - 2) confounding variables may be specified.}
  \item{k}{The fold of the cross validation. If k is NULL then only the first stage of MOSS is carried out.}  
}
\value{A list with 4 data frame elements:
\item{topRegressions}{The top regressions identified together with their log marginal likelihood}
\item{postIncProbs}{The posterior inclusion probabilities of each SNP that appears in one of the top regressions. This is obtained by adding the marginal likelihoods of the regressions in which each SNP appears and then normalizing over all the regressions.}
\item{interactionModels}{The best (in terms of marginal likelihood) hierarchical log-linear model containing the variables in each of the top regressions.}
\item{fits}{The fitted interaction models (using the glm function).}
\item{cvMatrix}{A matrix with the overall results of the k-fold cross validation. This table is typically called a confusion matrix.}
\item{cvDiag}{Some diagnostic information based on the cross validation: 'acc' is the accuracy, 'tpr' is the true positive rate, 'fpr' is the false postive rate, and 'auc' is the area under the ROC curve.}
}
\references{[1] Massam, H., Liu, J. and Dobra, A. (2009). A conjugate prior for discrete hierarchical log-linear models. Annals of Statistics, 37, 3431-3467.

[2] Dobra, A., Briollais, L., Jarjanazi, H., Ozcelik, H. and Massam, H. (2010). Applications of the mode oriented stochastic search (MOSS) algorithm for discrete multi-way data to genomewide studies. Bayesian Modeling in Bioinformatics, Taylor & Francis (Dey, D., Ghosh, S., and Mallick, B., eds.), 63-93. 

[3] Dobra, A. and Massam, H. (2010). The mode oriented stochastic search (MOSS) algorithm for log-linear models with conjugate priors. Statistical Methodology, 7, 240-253. 
}
\author{Matthew Friedlander, Adrian Dobra, Helene Massam, and Laurent Briollais} 
\note{The function recode_data can decide whether diallelic SNPs (i.e. SNPs with three categories) should be recoded as binary and in which way. If desired this function should be used prior to MOSS_GWAS.}
\seealso{\code{\link{recode_data}}}
\examples{
data(simuCC)
data <- simuCC[,c(1002,2971,rep(5978:6001))]
# The SNPs in columns 1002 and 2971 of simuCC called rs4491689 and rs6869003 cause the disease.
set.seed(7)
MOSS_GWAS (alpha = 1, c = 0.1, cPrime = 0.0001, q = 0.1, replicates = 1, 
           maxVars = 3, data, dimens = c(rep(3,25),2), confVars = NULL, k = NULL)
}
\keyword{htest}
\keyword{models}
