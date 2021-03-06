\name{EstimateTFDR}
\alias{EstimateTFDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate T False Discovery Rates.}
\description{
An estimated false discovery rate for each gene is computed, based on
the T scores. The T scores are as follows, for two-class, survival, and
quantitative outcomes: two-sample t-statistics, Cox scores, standardized
regression coefficients. The output of this function is identical to the
outputs "fdrt" and "pi0" of the function EstimateLPCFDR. This function
should be used if only the FDR of T is desired, because computing the
FDR of LPC is time-consuming.
}
\usage{
EstimateTFDR(x,y, type,censoring.status=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The matrix of gene expression values; pxn where n is the
    number of observations and p is the number of genes.}
  \item{y}{A vector of length n, with an outcome for each observation. For two-class
   outcome, y's elements are 1 or 2. For quantitative outcome, y's
   elements are real-valued. For survival
  data, y indicates the survival time. For a multiclass outcome, y is
  coded as 1,2,3,...}    
  \item{type}{One of "regression" (for a quantitative outcome), "two
    class", "survival", or "multiclass" (for a multiple-class outcome).}
  \item{censoring.status}{For survival outcome only, a vector of
      length
      n which takes on values 0 or 1
      depending on whether the observation is complete or censored.}
}
\details{
False discovery rates are estimated by permutations, as in e.g. Tusher
et al (2001) and Storey & Tibshirani (2003).
}
\value{
 % ~Describe the value returned
 % If it is a LIST, use
  \item{fdrt}{A vector of length p (equal to the number of genes),
    with the T false discovery rate for each gene. Note that this is
    identical to the "fdrt" output by the function EstimateLPCFDR.} 
  \item{pi0}{The fraction of genes that are believed to be null.}
  \item{call}{The function call made.}
}
\references{
Storey, J.D. and Tibshirani, R. (2003) Statistical significance for
genomewide studies. Proceedings of the National Academy of
Sciences. 100(16): 9440-9445.

Tusher, V.G. and Tibshirani, R. and Chu, G. (2001) Significance analysis
of microarrays applied to the ionizing radiation response. Proceedings
of the National Academy of Sciences. 98(9): 5116-5121.

  Witten, D.M. and Tibshirani, R. (2008) Testing significance
  of features by lassoed principal components. Annals of Applied
  Statistics. http://www-stat.stanford.edu/~dwitten}
\author{Daniela M. Witten and Robert Tibshirani}
\examples{
### not running due to timing - uncomment to run

#set.seed(2)
#n <- 40 # 40 samples
#p <- 1000 # 1000 genes
#x <- matrix(rnorm(n*p), nrow=p) # make 40x1000 gene expression matrix
#y <-  rnorm(n) # quantitative outcome
## make first 50 genes differentially-expressed
#x[1:25,y<(-.5)] <- x[1:25,y<(-.5)]+ 1.5
#x[26:50,y<0] <- x[26:50,y<0] - 1.5
## compute LPC and T scores for each gene
#lpc.obj <- LPC(x,y, type="regression")
## Look at plot of Predictive Advantage
#pred.adv <-
#PredictiveAdvantage(x,y,type="regression",soft.thresh=lpc.obj$soft.thresh)
## Estimate FDRs for LPC and T scores
#fdr.lpc.out <-
#EstimateLPCFDR(x,y,type="regression",soft.thresh=lpc.obj$soft.thresh,nreps=50)
## Estimate FDRs for T scores only. This is quicker than computing FDRs
##    for LPC scores, and should be used when only T FDRs are needed. If we
##    started with the same random seed, then EstimateTFDR and EstimateLPCFDR
##    would give same T FDRs.
#fdr.t.out <- EstimateTFDR(x,y, type="regression")
## print out results of main function
#lpc.obj
## print out info about T FDRs
#fdr.t.out
## print out info about LPC FDRs
#fdr.lpc.out
## Compare FDRs for T and LPC on 6\% of genes. In this example, LPC has
##    lower FDR.
#PlotFDRs(fdr.lpc.out,frac=.06)
## Print out names of 20 genes with highest LPC scores, along with their
##    LPC and T scores.
#PrintGeneList(lpc.obj,numGenes=20)
## Print out names of 20 genes with highest LPC scores, along with their
##    LPC and T scores and their FDRs for LPC and T.
#PrintGeneList(lpc.obj,numGenes=20,lpcfdr.out=fdr.lpc.out)

# Now, repeating everything that we did before, but using a
#   **survival** outcome
# NOT RUNNING DUE TO TIMING -- UNCOMMENT TO RUN

#set.seed(2)
#n <- 40 # 40 samples
#p <- 1000 # 1000 genes
#x <- matrix(rnorm(n*p), nrow=p) # make 40x1000 gene expression matrix
#y <-  rnorm(n) + 10 # survival times; must be positive
## censoring outcome: 0 or 1
#cens <- rep(1,40) # Assume all observations are complete
## make first 50 genes differentially-expressed
#x[1:25,y<9.5] <- x[1:25,y<9.5]+ 1.5
#x[26:50,y<10] <- x[26:50,y<10] - 1.5
#lpc.obj <- LPC(x,y, type="survival", censoring.status=cens)
## Look at plot of Predictive Advantage
#pred.adv <-
#PredictiveAdvantage(x,y,type="survival",soft.thresh=lpc.obj$soft.thresh,
#censoring.status=cens)
## Estimate FDRs for LPC scores and T scores
#fdr.lpc.out <- EstimateLPCFDR(x,y,
#type="survival",soft.thresh=lpc.obj$soft.thresh,
#nreps=20,censoring.status=cens)
## Estimate FDRs for T scores only. This is quicker than computing FDRs
##    for LPC scores, and should be used when only T FDRs are needed. If we
##    started with the same random seed, then EstimateTFDR and EstimateLPCFDR
##    would give same T FDRs.
#fdr.t.out <- EstimateTFDR(x,y, type="survival", censoring.status=cens)
## print out results of main function
#lpc.obj
## print out info about T FDRs
#fdr.t.out
## print out info about LPC FDRs
#fdr.lpc.out
## Compare FDRs for T and LPC scores on 10\% of genes.
#PlotFDRs(fdr.lpc.out,frac=.1)
## Print out names of 20 genes with highest LPC scores, along with their
##    LPC and T scores.
#PrintGeneList(lpc.obj,numGenes=20)
## Print out names of 20 genes with highest LPC scores, along with their
##    LPC and T scores and their FDRs for LPC and T.
#PrintGeneList(lpc.obj,numGenes=20,lpcfdr.out=fdr.lpc.out)

}
