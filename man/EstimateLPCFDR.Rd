\name{EstimateLPCFDR}
\alias{EstimateLPCFDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate LPC False Discovery Rates.}
\description{
An estimated false discovery rate for each gene is computed, based on
the LPC scores. Estimated false discovery rates based on T scores (Cox,
regression coefficient, or two-sample t-statistic) are also given.
}
\usage{
EstimateLPCFDR(x, y, type, nreps=100,soft.thresh=NULL, censoring.status=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The matrix of gene expression values; pxn where n is the
    number of observations and p is the number of genes.}
  \item{y}{A vector of
    length n, with an outcome for each observation. For two-class outcome, y's
    elements are 1 or 2. For quantitative outcome, y's elements are
    real-valued. For survival
    data, y indicates the survival time. For multiclass, y is coded as 1,2,3,...}
  \item{type}{One of "regression" (for a quantitative outcome), "two
    class", "multiclass", or "survival".}
  \item{nreps}{The number of training/test set splits used to estimate
    LPC's false discovery rates. Default is 100.}
  \item{soft.thresh}{The value of the soft threshold to be used in the
    L1 regression of the scores onto the eigenarrays. This can be the
    value of the soft-threshold chosen adaptively by LPC function when it
    is run on the data; if not entered by the user, then that adaptive
    value is computed.}
  \item{censoring.status}{For survival outcome only, a vector of length
    n which takes on values 0 or 1
  depending on whether the observation is complete or censored.}
}
\details{
Details of false discovery rate estimation for LPC can be found in the
paper: http://www-stat.stanford.edu/~dwitten \\

As explained in the paper, FDR of LPC is estimated by computing the FDR
of simpler scores (Cox for survival outcome, standardized regression
coefficients for regression outcome, and two-sample t-statistic for
two-class outcome) and then estimating the difference between the FDR of LPC
and the FDR of these simpler scores. 
}
\value{
 % ~Describe the value returned
 % If it is a LIST, use
  \item{fdrlpc}{A vector of length p (equal to the number of genes),
    with the LPC false discovery rate for each gene.}
  \item{fdrt}{A vector of length p (equal to the number of genes),
    with the T false discovery rate for each gene.}
  \item{fdrdiff}{The FDR of T minus the FDR of LPC. This is
    approximately equal to fdrt-fdrlpc, with the caveat that the FDR of
    LPC (computed via this difference) must be between 0 and 1.}
  \item{pi0}{The fraction of genes that are believed to be null.}
  \item{soft.thresh}{The value of the soft threshold used in the L1
    constraint for LPC.}
  \item{call}{The function call made.}
}
\references{Witten, D.M. and Tibshirani, R. (2008) Testing significance
  of features by lassoed principal components. Annals of Applied
  Statistics. http://www-stat.stanford.edu/~dwitten}
\author{Daniela M. Witten and Robert Tibshirani}
\note{
}
\seealso{}
\examples{
set.seed(2)
n <- 40 # 40 samples
p <- 1000 # 1000 genes
x <- matrix(rnorm(n*p), nrow=p) # make 40x1000 gene expression matrix
y <-  rnorm(n) # quantitative outcome
# make first 50 genes differentially-expressed
x[1:25,y<(-.5)] <- x[1:25,y<(-.5)]+ 1.5
x[26:50,y<0] <- x[26:50,y<0] - 1.5
# compute LPC and T scores for each gene
lpc.obj <- LPC(x,y, type="regression")
# Look at plot of Predictive Advantage
pred.adv <-
PredictiveAdvantage(x,y,type="regression",soft.thresh=lpc.obj$soft.thresh)
# Estimate FDRs for LPC and T scores
fdr.lpc.out <-
EstimateLPCFDR(x,y,type="regression",soft.thresh=lpc.obj$soft.thresh,nreps=50)
# Estimate FDRs for T scores only. This is quicker than computing FDRs
#    for LPC scores, and should be used when only T FDRs are needed. If we
#    started with the same random seed, then EstimateTFDR and EstimateLPCFDR
#    would give same T FDRs.
fdr.t.out <- EstimateTFDR(x,y, type="regression")
# print out results of main function
lpc.obj
# print out info about T FDRs
fdr.t.out
# print out info about LPC FDRs
fdr.lpc.out
# Compare FDRs for T and LPC on 6\% of genes. In this example, LPC has
#    lower FDR.
PlotFDRs(fdr.lpc.out,frac=.06)
# Print out names of 20 genes with highest LPC scores, along with their
#    LPC and T scores.
PrintGeneList(lpc.obj,numGenes=20)
# Print out names of 20 genes with highest LPC scores, along with their
#    LPC and T scores and their FDRs for LPC and T.
PrintGeneList(lpc.obj,numGenes=20,lpcfdr.out=fdr.lpc.out)

# Now, repeating everything that we did before, but using a
#   **survival** outcome

set.seed(2)
n <- 40 # 40 samples
p <- 1000 # 1000 genes
x <- matrix(rnorm(n*p), nrow=p) # make 40x1000 gene expression matrix
y <-  rnorm(n) + 10 # survival times; must be positive
# censoring outcome: 0 or 1
cens <- rep(1,40) # Assume all observations are complete
# make first 50 genes differentially-expressed
x[1:25,y<9.5] <- x[1:25,y<9.5]+ 1.5
x[26:50,y<10] <- x[26:50,y<10] - 1.5
lpc.obj <- LPC(x,y, type="survival", censoring.status=cens)
# Look at plot of Predictive Advantage
pred.adv <- PredictiveAdvantage(x,y,type="survival",soft.thresh=lpc.obj$soft.thresh,
censoring.status=cens)
# Estimate FDRs for LPC scores and T scores
fdr.lpc.out <- EstimateLPCFDR(x,y,type="survival",
soft.thresh=lpc.obj$soft.thresh,nreps=20,censoring.status=cens)
# Estimate FDRs for T scores only. This is quicker than computing FDRs
#    for LPC scores, and should be used when only T FDRs are needed. If we
#    started with the same random seed, then EstimateTFDR and EstimateLPCFDR
#    would give same T FDRs.
fdr.t.out <- EstimateTFDR(x,y, type="survival", censoring.status=cens)
# print out results of main function
lpc.obj
# print out info about T FDRs
fdr.t.out
# print out info about LPC FDRs
fdr.lpc.out
# Compare FDRs for T and LPC scores on 10\% of genes.
PlotFDRs(fdr.lpc.out,frac=.1)
# Print out names of 20 genes with highest LPC scores, along with their
#    LPC and T scores.
PrintGeneList(lpc.obj,numGenes=20)
# Print out names of 20 genes with highest LPC scores, along with their
#    LPC and T scores and their FDRs for LPC and T.
PrintGeneList(lpc.obj,numGenes=20,lpcfdr.out=fdr.lpc.out)

}
\keyword{}
\keyword{}% __ONLY ONE__ keyword per line
