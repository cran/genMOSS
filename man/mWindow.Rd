\name{mWindow}
\Rdversion{1.2}
\title{Analyzing GWAS data in sequence over a moving window}
\description{A simple alternative to stochastically searching through all combinations of SNPs (as MOSS_GWAS does) is to group SNPs together according to the sequence that they appear in a genetic region (if that information is available). Definining a window size w we first group SNPs 1 to w together and then SNPs 2 to w + 1 together and so on. If we let Y be a response, X be a group of SNPs, and r = Y|X be called the regression of Y on X, the mWindow function computes the marginal likelihood, P(r) = P(Y|X), of the regression corresponding to each group (or window). The aim is to identify those regressions such that the marginal likelihood is the highest. The groups of SNPs (or genetic regions) contained in these regressions are most associated with the response. The prior used in Bayesian computations is the generalized hyper Dirichlet of Massam et. al (2009).}
\usage{mWindow (data, dimens, alpha = 1, windowSize = 2)}
\alias{mWindow}
\arguments{
  \item{data}{A data frame containing the genotype information for a given set of SNPs. The data frame should be organized such that each row refers to a subject and each column to a SNP. The last column must be a binary response for each subject. The data frame must contain at least 8 columns. Rows containing any missing values (i.e. NAs) are omitted from the analysis.} 
  \item{dimens}{The number of possible values for each column of data. Each possible value does not need to occur in data. All entries of dimens must be greater than or equal to 2.} 
  \item{alpha}{A hyperparameter of the prior representing the total of a fictive contingency table with counts equal to alpha divided by the number of cells. Alpha must be a positive real number.}
  \item{windowSize}{The size of the moving window. Must be an integer from 1 to 5.}
}
\value{A data frame listing the regression formed in each window and its corresponding log marginal likelihood. The data frame is sorted in descending order by log marginal likelihood.}
\references{[1] Massam, H., Liu, J. and Dobra, A. (2009). A conjugate prior for discrete hierarchical log-linear models. Annals of Statistics, 37, 3431-3467.

[2] Dobra, A., Briollais, L., Jarjanazi, H., Ozcelik, H. and Massam, H. (2010). Applications of the mode oriented stochastic search (MOSS) algorithm for discrete multi-way data to genomewide studies. Bayesian Modeling in Bioinformatics, Taylor & Francis (Dey, D., Ghosh, S., and Mallick, B., eds.), 63-93. 

[3] Dobra, A. and Massam, H. (2010). The mode oriented stochastic search (MOSS) algorithm for log-linear models with conjugate priors. Statistical Methodology, 7, 240-253. 

[4] Sun, J., Levin, A., Boerwinkle, E., Robertson, H., and Kardia, S. (2006). A scan statistic for identifying chromosomal patterns of SNP association. Genetic Epidemiology, 30, 627-635.

[5] Wu, M., Kraft, P., Epstein, M., Taylor, D., Chanock, S., Hunter, D., and Lin, X. (2010). Powerful SNP-set analysis for case-control genome-wide association studies. The American Journal of Human Genetics, 86, 929-942.
}
\author{Matthew Friedlander, Adrian Dobra, Helene Massam, and Laurent Briollais} 
\note{The function recode_data can decide whether diallelic SNPs (i.e. SNPs with three categories) should be recoded as binary and in which way. If desired this function should be used prior to mWindow.}
\seealso{\code{\link{recode_data}}}
\examples{
data(simuCC)
data <- simuCC[,c(1002,2971,rep(5978:6001))]
# The SNPs in columns 1002 and 2971 of simuCC called rs4491689 and rs6869003 cause the disease.
set.seed(7)
s <- mWindow (data, dimens = c(rep(3,25),2), alpha = 1, windowSize = 2)
head (s, n = 5)
}
\keyword{htest}
\keyword{models}
