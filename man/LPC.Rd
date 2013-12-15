\name{LPC}
\alias{LPC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute LPC score for each gene}
\description{
This is the main function for the Lassoed Principal Components (lpc)
method. Given a matrix of gene expression data and a vector of outcome
  measurements (one per patient) and an outcome type, this function
  computes the LPC score for each gene. Simple scores
  are projected onto the eigenarrays (with an L1 constraint) in
  order to obtain LPC scores; these simple scores are Cox scores (for a survival outcome),
  two-sample t-statistics (for a two-class outcome) or standardized
  regression coefficients (for a regression/quantitative
  outcome).
}
\usage{
LPC(x, y, type = "regression", soft.thresh = NULL, u =
NULL, censoring.status = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The matrix of gene expression values; pxn where n is the
    number of observations and p is the number of genes.}
  \item{y}{A vector of
        length n, with an outcome for each observation. For two-class
        outcome, y's
            elements are 1 or 2. For quantitative outcome, y's elements
            are
            real-valued. For survival
                data, y indicates the survival time, and must be
                positive. For multiclass outcome, y must be coded as 1,2,3,...}   
  \item{type}{One of "regression" (for a quantitative outcome), "two
    class", "multiclass", or "survival".}
  \item{soft.thresh}{The value of the soft threshold to be used in the
    L1 regression of the scores onto the eigenarrays. If NULL, then it
    will be chosen adaptively.}
  \item{u}{The eigenarrays of the data matrix (useful if the data matrix
  is large, to avoid having to re-compute the SVD in order to run this
  function multiple times). If NULL, then the eigenarrays will be
  computed internally. Note that the eigenarrays have the same dimension
  as the data.}
\item{censoring.status}{For survival outcome only, a vector of length
  n which takes on values 0 or 1
    depending on whether the observation is complete or censored.}      
}
\details{
Lassoed Principal Components (LPC) is a method for identifying
features associated with an outcome, under the hypothesis that
significant features occur in correlated sets. In particular, the method
was designed for identifying genes in a microarray experiment that are
associated with some outcome (e.g.: cancer vs. normal, survival time,
tumor size). The method involves projecting simpler scores for each gene
(e.g.: two-sample t-statistic, Cox score, standardized regression
coefficients) onto the eigenarrays of the gene expression data matrix,
subject to an L1 constraint. The fitted values are the LPC scores, which
may result in a more accurate gene ranking than the initial simpler
scores. 
}
\value{
 % ~Describe the value returned
 % If it is a LIST, use
 \item{lpcscores}{A vector of LPC scores, of length p (one for each
   gene).}
 \item{tscores}{A vector of T scores (Cox, two-sample t-stats, or
   standardized regression coefficients), of length p (one for each
   gene). }
  \item{soft.thresh}{The value of the soft threshold used.}
  \item{coefs}{The coefficients for the L1 regression of the scores onto
    the eigenarrays.}
  \item{call}{The call that was made to this function.}
}
\references{Witten, D.M. and Tibshirani, R. (2008) Testing significance
  of features by lassoed principal components. Annals of Applied
  Statistics. http://www-stat.stanford.edu/~dwitten}
\author{Daniela M. Witten and Robert Tibshirani}
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
pred.adv <- PredictiveAdvantage(x,y,type="regression",soft.thresh=lpc.obj$soft.thresh)
# Estimate FDRs for LPC and T scores
fdr.lpc.out <-
EstimateLPCFDR(x,y,type="regression",soft.thresh=lpc.obj$soft.thresh,nreps=5)
# should run more reps in practice! Running just 5 in this example
# so that it runs fast.
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
### COMMENTED TO REDUCE RUN TIME -- UNCOMMENT TO RUN

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
#fdr.lpc.out <-
#EstimateLPCFDR(x,y, type="survival",soft.thresh=lpc.obj$soft.thresh,
#nreps=5,censoring.status=cens)
# Running just 5 reps so that things run quickly -- please run more in practice
## Estimate FDRs for T scores only. This is quicker than computing FDRs
###    for LPC scores, and should be used when only T FDRs are needed. If we
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


## Now, repeating everything that we did before, but using a
##   **multiclass** outcome

#set.seed(2)
#n <- 40 # 40 samples
#p <- 1000 # 1000 genes
#x <- matrix(rnorm(n*p), nrow=p) # make 40x1000 gene expression matrix
#y <-  sample(1:4, n, rep=TRUE)
## make first 50 genes differentially-expressed
#x[1:50,y==1] <- x[1:50,y==1] + 1.5
#x[1:25,y==2] <- x[1:25,y==2] + 1
#x[1:50,y==3] <- x[1:50,y==3] - 1
#lpc.obj <- LPC(x,y, type="multiclass")
## Look at plot of Predictive Advantage
#pred.adv <-
#PredictiveAdvantage(x,y,type="multiclass",soft.thresh=lpc.obj$soft.thresh)
## Estimate FDRs for LPC scores and T scores
#fdr.lpc.out <-
#EstimateLPCFDR(x,y, type="multiclass",soft.thresh=lpc.obj$soft.thresh,
#nreps=20)
## Estimate FDRs for T scores only. This is quicker than computing FDRs
##    for LPC scores, and should be used when only T FDRs are needed. If we
##    started with the same random seed, then EstimateTFDR and EstimateLPCFDR
##    would give same T FDRs.
#fdr.t.out <- EstimateTFDR(x,y, type="multiclass")
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
