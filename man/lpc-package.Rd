\name{lpc-package}
\alias{lpc-package}
\alias{lpc}
\docType{package}
\title{
  This package implements the Lassoed Principal Components (LPC) method.
  It is used to compute LPC scores for each gene in a microrray
  experiment with a survival, quantitative, or two-class outcome.
}
\description{
  LPC scores are computed for each gene; the method borrows strength
  across genes and can result in more accurate gene scores than simpler
  statistics. In this package, the LPC method is applied by regressing Cox scores
  (survival outcome), two-sample t-statistics (two-class outcome), or
  standardized regression coefficients (quantitative outcome) onto gene
  expression eigenarrays, with an L1 constraint.
}
\details{
\tabular{ll}{
Package: \tab lpc\cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2013-12-15\cr
License: \tab GPL (>= 2)
}
The main function is "LPC", which computes LPC scores for each
gene. The matrix of gene expression data, a vector of outcome, and the
outcome type must be passed in to this function.
}
\author{
  Daniela M. Witten and Robert Tibshirani
  
Maintainer: Daniela M. Witten <dwitten@uw.edu>
}
\references{
  Witten, DM and Tibshirani R (2008) Testing significance of features by
  lassoed principal components. Annals of Applied Statistics.
}

\keyword{ package }
\seealso{
www.biostat.washington.edu/~dwitten
}
\examples{
### Uncomment to run


#n <- 40 # 40 samples
#p <- 1000 # 1000 genes
#x <- matrix(rnorm(n*p), nrow=p) # make 40x1000 gene expression matrix
#y <-  rnorm(n) # quantitative outcome
## make first 50 genes differentially-expressed
#x[1:25,y<(-.5)] <- x[1:25,y<(-.5)]+ 1.5
#x[26:50,y<0] <- x[26:50,y<0] - 1.5
## compute LPC and T scores for each gene
#lpc.obj <- LPC(x,y, type="regression")
#lpc.obj
## Look at plot of Predictive Advantage
#pred.adv <- PredictiveAdvantage(x,y,type="regression",soft.thresh=lpc.obj$soft.thresh)
## Estimate FDRs for LPC and T scores
#fdr.lpc.out <- EstimateLPCFDR(x,y,type="regression",soft.thresh=lpc.obj$soft.thresh,nreps=50)
#fdr.lpc.out
}
