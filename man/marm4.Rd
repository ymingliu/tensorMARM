\name{marm4}
\alias{marm4-function}
\alias{marm4}
\docType{package}
\title{
  Fit structural MARM with sparsity assumption and fixed ranks.
}
\description{
  Fit a structural intergative multi-view multivariate additive model (structural marm) using B-splines with given ranks (\eqn{r_{1}, r_{2}, r_{3}, r_{4}}). A fourth-order coefficient tensor can be estimated by this function. The group sparse penalty such as \code{LASSO}, \code{MCP} or \code{SCAD} and the coordinate descent algorithm are used to yield a sparsity estimator. The \code{BIC} or \code{cross-validation} method are used to search the optimal regularization parameter.
}

\usage{
marm4 <- function(Y,X,group=NULL,K=6,r1=NULL,r2=NULL,r3=NULL,r4=NULL,
                   method="BIC",ncv=10,penalty="LASSO",lambda=NULL,D0=NULL,
                   intercept=TRUE,nlam=20,degr=3,lam_min=0.01, 
                   eps=1e-4,max_step=10,eps1=1e-4,max_step1=10,
                   gamma=2,dfmax=NULL,alpha=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times p} numeric design matrix for the model, where \eqn{p=\sum_{g}p_g}}
  
  \item{group}{A \eqn{p} vector of the grouping index of predictors, e.g., \eqn{group=c(1,1,1,2,2,2)} means there are \eqn{6} perdictors in the model, and the first three perdictors are in the same group and the last three perdictors are in another one. By default, we set \eqn{group=rep(1,p)}.}
  
  \item{K}{The number of B-spline basis functions, that is the plus of both degrees of basis functioins and the number of knots. Default is \code{6}, which means cubic splines.}
  
  \item{r1}{The first dimension of single value matrix of the tensor. Default is \code{2}.}
  
  \item{r2}{The second dimension of single value matrix of the tensor. Default is \code{2}.}
  
  \item{r3}{The third dimension of single value matrix of the tensor. Default is \code{2}.}
  
  \item{r4}{The fourth dimension of single value matrix of the tensor. Default is \code{2}.}
  
  \item{method}{The method to be applied to select regularization parameters.  Either \code{BIC} (default), or \code{CV}.}
    
  \item{ncv}{The number of cross-validation folds. Default is \code{10}. If \code{method} is not \code{CV}, \code{ncv} is useless.}  
  
  \item{penalty}{The penalty to be applied to the model. Either \code{LASSO} (the default), \code{MCP} or \code{SCAD}.}
  
  \item{lambda}{A user-specified sequence of lambda values.  By default, a sequence of values of length \code{nlam} is computed, equally spaced on the log scale.}
  
  \item{D0}{A user-specified list of initialized values, including five initialized matrix \eqn{S_{(4)}}  (called \code{S}), \code{A}, \code{B}, \code{C} and \code{D}. By default, a list of initialization satisfying fixed ranks is computed by random.}
 
  \item{intercept}{A logical value indicating whether the intercept is fitted. Default is \code{TRUE} or set to zero by \code{FALSE}.}
  
  \item{degr}{The number of knots of B-spline base function. Default is \code{3}.}
  
  \item{nlam}{The number of lambda values. Default is \code{20}.}
  
  \item{lam_min}{The smallest value for lambda, as a fraction of lambda.max.  Default is \code{0.01}.}
  
  \item{eps}{Convergence threshhold. The algorithm iterates until the relative change in any coefficient is less than \code{eps}. Default is \code{1e-4}.}
  
  \item{max_step}{Maximum number of iterations. Default is \code{20}.}
              
  \item{eps1}{Convergence threshhold. The Coordinate descent method algorithm iterates until the relative change in any coefficient is less than \code{eps1}. Default is \code{1e-4}.}
  
  \item{max_step1}{The maximum iterates number of coordinate descent method. Default is \code{20}.}
  
  \item{gamma}{The tuning parameter of the MCP/SCAD penalty.}
  
  \item{dfmax}{Upper bound for the number of nonzero coefficients. Default is no upper bound.  However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.}
               
  \item{alpha}{Tuning parameter for the Mnet estimator which controls the relative contributions from the LASSO, MCP/SCAD penalty and the ridge, or L2 penalty.  \code{alpha=1} is equivalent to LASSO, MCP/SCAD penalty, while \code{alpha=0} would be equivalent to ridge regression. However, \code{alpha=0} is not supported; \code{alpha} may be arbitrarily small, but not exactly 0.}
  
}

\details{
  This function gives \code{pq} functional coefficients' estimators of marm. A fourth-order tensor with multiple ranks (\eqn{r_1, r_2, r_3, r_4}) need to be estimated. We fix these ranks and use an alternative updating algorithm to update its core tensor and factor matrices based on Tucker decomposition. Group LASSO, SCAD or MCP penalty is applied on the row of the factor matrix \eqn{A} to achieve variable selection.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{D}{Estimator of \eqn{D_{(4)}}.}
  
  \item{mu}{Estimator of intercept \eqn{\mu}.}
  
  \item{S.opt}{A length-\eqn{ng} list including estimator of the core tensor \eqn{S_{(3)}} of each coefficient tensor.}
  
  \item{A.opt}{A length-\eqn{ng} list including estimator of the factor matrix \eqn{A} of each coefficient tensor.}
  
  \item{B.opt}{A length-\eqn{ng} list including estimator of the factor matrix \eqn{B} of each coefficient tensor.}
  
  \item{C.opt}{A length-\eqn{ng} list including estimator of the factor matrix \eqn{C} of each coefficient tensor.}
  
  \item{D.opt}{A length-\eqn{ng} list including estimator of the factor matrix \eqn{D} of each coefficient tensor.}
  
  \item{lambda.seq}{The sequence of regularization parameter values in the path.}
  
  \item{lambda_opt}{The value of \code{lambda} with the minimum \code{BIC} or \code{CV} value.}
  
  \item{rss}{Residual sum of squares (RSS).}
  
  \item{df}{Degrees of freedom.}
  
  \item{activeX}{The active set of \eqn{X}. A length-\eqn{p} vector.}
  
  \item{opts}{Other related parameters used in algorithm. Some of them are set by default.}
  
  \item{opts}{Other related parameters used in algorithm (especially parameters in penalty). Some of them are set by default.}

  %\item{...}{ Other options for CompositeQuantile.}
}

%\author{
%Your Name, email optional.
%Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
%}
\references{
  Multivariate additive regression for multi-view data via tensor estimation.
}
\keyword{ Group sparsity; Tensor low-rankness; Tucker decomposition; Multivariate additive regression; Multi-view data. }
\seealso{
  marm4_dr
}

\examples{
library(tensorMARM)
n <- 200; q <- 5; p <- 100; s <- 3; ng <- 4
group <- rep(1:ng,each=p/ng)
mydata <- marm4.sim.fbs(n,q,p,s,group,isfixedR=1)
fit <- with(mydata, marm4(Y,X,group,K,r10,r20,r30,r40,D0=D0,nlam=5))
}