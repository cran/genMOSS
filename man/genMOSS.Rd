\name{genMOSS}
\Rdversion{1.2}
\title{Functions for the Bayesian analysis of GWAS data}
\description{The genMOSS package implements the Mode Oriented Stochastic Search (MOSS) of Dobra and Massam (2010) as well as a simple moving window approach (see Sun et. al (2006) and Wu et. al (2010)) to identify combinations of SNPs that are associated with a response.} 
\details{
\tabular{ll}{
Package: \tab genMOSS\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2014-12-01\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The function MOSS_GWAS implements the MOSS algorithm while the mWindow implements the moving window approach.}
\author{Authors: Matthew Friedlander, Adrian Dobra, Helene Massam, and Laurent Briollais,
Maintainer: Matthew Friedlander <friedlander.matthew@gmail.com>
}
\references{[1] Massam, H., Liu, J. and Dobra, A. (2009). A conjugate prior for discrete hierarchical log-linear models. Annals of Statistics, 37, 3431-3467.

[2] Dobra, A., Briollais, L., Jarjanazi, H., Ozcelik, H. and Massam, H. (2010). Applications of the mode oriented stochastic search (MOSS) algorithm for discrete multi-way data to genomewide studies. Bayesian Modeling in Bioinformatics, Taylor & Francis (Dey, D., Ghosh, S., and Mallick, B., eds.), 63-93.

[3] Dobra, A. and Massam, H. (2010). The mode oriented stochastic search (MOSS) algorithm for log-linear models with conjugate priors. Statistical Methodology, 7, 240-253.

[4] Sun, J., Levin, A., Boerwinkle, E., Robertson, H., and Kardia, S. (2006). A scan statistic for identifying chromosomal patterns of SNP association. Genetic Epidemiology, 30, 627-635.

[5] Wu, M., Kraft, P., Epstein, M., Taylor, D., Chanock, S., Hunter, D., and Lin, X. (2010). Powerful SNP-set analysis for case-control genome-wide association studies. The American Journal of Human Genetics, 86, 929-942. 
}
\examples{
data(simuCC)
data <- simuCC[,c(1002,2971,rep(5978:6001))]
# The SNPs in columns 1002 and 2971 of simuCC called rs4491689 and rs6869003 cause the disease.
set.seed(7)
MOSS_GWAS (alpha = 1, c = 0.1, cPrime = 0.0001, q = 0.1, replicates = 1, 
           maxVars = 3, data, dimens = c(rep(3,25),2), confVars = NULL, k = NULL)
s <- mWindow (data, dimens = c(rep(3,25),2), alpha = 1, windowSize = 2)
head (s, n = 5)
}
\keyword{htest}
\keyword{models}
