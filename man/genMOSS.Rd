\name{genMOSS}
\Rdversion{1.0}
\title{An implementation of the MOSS algorithm for the analysis of GWAS data.}
\description{The MOSS algorithm is a Bayesian variable selection procedure that can be used for the analysis GWAS data. It identifies combinations of the best predictive SNPs associated with the response. It also performs a hierarchical log-linear model search to identify the most relevant associations among the resulting subsets of SNPs. The prior used is the generalized hyper Dirichlet.}
\details{
\tabular{ll}{
Package: \tab genMOSS\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-03-06\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The package contains the function MOSS.GWAS which carries out the MOSS search. The function has an option to use model averaging to construct a classifier for predicting the response and to assess its capability via k-fold cross validation.}
\author{Authors: Matthew Friedlander and Laurent Briollais, 
Maintainer: Matthew Friedlander <friedlander@lunenfeld.ca>
}
\references{[1] Massam, H., Liu, J. and Dobra, A. (2009). A conjugate prior for discrete hierarchical log-linear models. Annals of Statistics, 37, 3431-3467.

[2] Dobra, A., Briollais, L., Jarjanazi, H., Ozcelik, H. and Massam, H. (2010). Applications of the mode oriented stochastic search (MOSS) algorithm for discrete multi-way data to genomewide studies. Bayesian Modeling in Bioinformatics, Taylor & Francis (D. Dey, S. Ghosh and B. Mallick, eds.), 63-93.

[3] Dobra, A. and Massam, H. (2010). The mode oriented stochastic search (MOSS) algorithm for log-linear models with conjugate priors. Statistical Methodology, 7, 240-253.
}
\examples{
data(simuCC)
data <- simuCC[,c(1002,2971,rep(5978:6001))]
# The SNPs in columns 1002 and 2971 of simuCC called rs4491689 and rs6869003 
# cause the disease.
set.seed(7)
MOSS.GWAS (alpha = 1, c = 0.1, cPrime = 0.0001, q = 0.1, replicates = 1, 
           maxVars = 3, data, dimens = c(rep(3,25),2), k = NULL)}
\keyword{htest}
\keyword{models}
