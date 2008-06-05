\name{PredictiveAdvantage}
\alias{PredictiveAdvantage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot predictive advantage of LPC vs. T}
\description{
This function plots the predictive advantage of LPC vs. T.  The
predictive advantage is the difference between the red and black
curves. If the red curve is higher than the black curve on average, then
LPC should be used instead of T on this data set. 
}
\usage{
PredictiveAdvantage(x,y,type,nreps=20,ngenes=100,soft.thresh=NULL,censoring.status=NULL)
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
            data, y indicates the survival time.}
  \item{type}{One of "regression" (for a quantitative outcome), "two
    class", or "survival".}
    \item{nreps}{Number of training/test set splits used in computing
      predictive advantage.}
    \item{soft.thresh}{Value of soft threshold used in L1 constraint for
      LPC. If NULL, then it will be computed adaptively.}
  \item{ngenes}{Number of genes to include in predictive advantage
    plot. (E.g., if ngenes=100 (default) then  $E(|T_test| |
    |L_train|>alpha(L_train))$ and 
$E(|T_test| | |T_train|>alpha(T_train))$ will be plotted (see
"Details"), where $alpha(L_train)$ and $alpha(T_train)$ are the 100th
largest (in absolute value) T and LPC scores on the training set.)}
\item{censoring.status}{For survival outcome only, a vector of length
  n which takes on values 0 or 1
    depending on whether the observation is complete or censored.}      
}
\details{
  As explained
in the paper, predictive advantage is computed by first splitting the
data into a training set and a test set (each with 50\% of the
samples). Then, the following is computed: $E(|T_test| | |L_train|>alpha(L_train)) -
E(|T_test| | |T_train|>alpha(T_train))$, where $T_test$ are the T scores on
the test data, $T_train$ are the T scores on the training data, and
$L_train$ are the LPC scores on the training data. $alpha(L_train)$ and
$alpha(T_train)$ are the alpha quantiles of the LPC and T scores on the training data. A
large value of $E(|T_test| | |L_train|>alpha(L_train)) -
E(|T_test| | |T_train|>alpha(T_train))$ suggests that LPC is superior to
T on this data set.

}
\value{
  \item{lpc}{A vector of numGenes elements. The ith element is of the form
    $E(|T_test| | |L_train|>alpha(L_train))$, where $alpha(L_train)$ is
    the ith smallest (in absolute value) training set LPC score.  Note that $L_train$
    are LPC scores on the training set, $T_test$ are T scores on the
    test set.}
  \item{t}{A vector of numGenes elements. The ith element is of the form
    $E(|T_test| | |L_train|>alpha(T_train))$, where $alpha(T_train)$ is
    the ith smallest (in absolute value) training set T score, and where $T_train$ are
  T scores on the training set and $T_test$ are T scores on the test set.}
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
pred.adv <-
PredictiveAdvantage(x,y,type="survival",soft.thresh=lpc.obj$soft.thresh,
censoring.status=cens)
# Estimate FDRs for LPC scores and T scores
fdr.lpc.out <-
EstimateLPCFDR(x,y,
type="survival",soft.thresh=lpc.obj$soft.thresh,nreps=20,censoring.status=cens)
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
