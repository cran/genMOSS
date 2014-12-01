\name{recode_data}
\Rdversion{1.2}
\title{Finding the optimal coding of diallelic SNPs}
\description{Let Y be a response, X be a SNP, and r = Y|X be called the regression of Y on X. For a diallelic SNP (i.e. a SNP with 3 categories), it may be that the marginal likelihood of the regression, P(r) = P(Y|X), is higher when the SNP is recoded as binary. Using the coding that maximizes this marginal likelihood may increase the power. Trinary variables can be recoded as binary in three different ways (or can be left as is). The function recode_data finds the optimal coding for each diallelic SNP in a given data frame and returns a revised data frame in the same order as the original. SNPs that are not diallelic are inserted into the new data frame unchanged. A vector containing the dimension of each SNP in the revised data frame is also returned. The prior used in Bayesian computations is the generalized hyper Dirichelt of Massam et. al (2009).} 
\usage{recode_data (data, dimens, alpha = 1)}
\alias{recode_data}
\arguments{
  \item{data}{A data frame containing the genotype information for a given set of SNPs. The data frame should be organized such that each row refers to a subject and each column to a SNP. The last column must be a binary response for each subject. The data frame must contain at least 8 columns. Rows containing any missing values (i.e. NAs) are omitted.} 
  \item{dimens}{The number of possible values for each column of data. Each possible value does not need to occur in data. All entries of dimens must be greater than or equal to 2.} 
  \item{alpha}{A hyperparameter of the prior representing the total of a fictive contingency table with counts equal to alpha divided by the number of cells. Alpha must be a positive real number.}
}

\value{A list with a data frame and a vector:
\item{recoded_data}{The recoded dataset.}
\item{recoded_dimens}{The revised dimension vector.}
}
\references{[1] Massam, H., Liu, J. and Dobra, A. (2009). A conjugate prior for discrete hierarchical log-linear models. Annals of Statistics, 37, 3431-3467.

[2] Dobra, A., Briollais, L., Jarjanazi, H., Ozcelik, H. and Massam, H. (2010). Applications of the mode oriented stochastic search (MOSS) algorithm for discrete multi-way data to genomewide studies. Bayesian Modeling in Bioinformatics, Taylor & Francis (Dey, D., Ghosh, S., and Mallick, B., eds.), 63-93. 

[3] Dobra, A. and Massam, H. (2010). The mode oriented stochastic search (MOSS) algorithm for log-linear models with conjugate priors. Statistical Methodology, 7, 240-253. 
}
\author{Matthew Friedlander, Adrian Dobra, Helene Massam, and Laurent Briollais} 
\examples{
data(simuCC)
data <- simuCC[,c(1002,2971,rep(5978:6001))]
# The SNPs in columns 1002 and 2971 of simuCC called rs4491689 and rs6869003 cause the disease.
set.seed(7)
r <- recode_data (data, dimens = c(rep(3,25),2), alpha = 1) 
s <- mWindow (data = r$recoded_data, dimens = r$recoded_dimens, alpha = 1, windowSize = 2)
head (s, n = 5)
}
\keyword{htest}
\keyword{models}
