\name{marm4.dr}
\alias{marm4.dr-function}
\alias{marm4.dr}
\docType{package}
\title{
  Fit structural MARM with sparsity assumption and unknown ranks.
}
\description{
  Fit a structural multivariate additive model for multi-view data (structural MARM) using B-splines with unknown ranks (\eqn{r_{1}, r_{2}, r_{3}, r_{4}}). A fourth-order coefficient tensor can be estimated by this function. The group sparse penalty such as LASSO, MCP or SCAD and the coordinate descent algorithm are used to yield a sparsity estimator. The BIC or cross-validation method is used to search the optimal regularization parameter, multiple ranks and the number of B-spline basis functions simultaneously.
}

\usage{
marm4.dr <- function(Y,X,group,K_index=NULL,r1_index=NULL,r2_index=NULL,
                      r3_index=NULL,r4_index=NULL,method="BIC",ncv=10,
                      penalty="LASSO",lambda=NULL,D0=NULL,
                      intercept=TRUE,nlam=20,degr=3,lam_min=0.01,
                      eps=1e-4,max_step=10,eps1=1e-4,max_step1=10,
                      gamma=2,dfmax=NULL,alpha=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{A \eqn{n\times q} numeric matrix of responses.}
  
  \item{X}{A \eqn{n\times p} numeric design matrix for the model, where \eqn{p=\sum_{g}p_g}}
  
  \item{group}{A \eqn{p} vector of the grouping index of predictors, e.g., \eqn{group=c(1,1,1,2,2,2)} means there are \eqn{6} perdictors in the model, and the first three perdictors are in the same group and the last three perdictors are in another one. By default, we set \eqn{group=rep(1,p)}.}
  
  \item{K_index}{The user-specified sequence of K. Default is a length-\code{1} vector \code{6}.}
  
  \item{r1_index}{A user-specified sequence of \eqn{r_1} values, where \eqn{r_1} is the first dimension of single value matrix of the tensor. Default is \code{r1_index}\eqn{=1,\cdots,\min(\lceil\log(n)\rceil,p)}.}
            
  \item{r2_index}{A user-specified sequence of \eqn{r_2} values, where \eqn{r_2} is the second dimension of single value matrix of the tensor. Default is \code{r2_index}\eqn{=1,\cdots},\code{max{K_index}}.}
                  
  \item{r3_index}{A user-specified sequence of \eqn{r_3} values, where \eqn{r_3} is the third dimension of single value matrix of the tensor. Default is \code{r3_index}\eqn{=1,\cdots,\min(\lceil\log(n)\rceil,ng)}.}
  
  \item{r4_index}{A user-specified sequence of \eqn{r_4} values, where \eqn{r_4} is the third dimension of single value matrix of the tensor. Default is \code{r4_index}\eqn{=1,\cdots,\min(\lceil\log(n)\rceil,q)}.}
  
  \item{method}{The method to be applied to search the number of B-spline basis functions, regularization parameters and multiple ranks simultaneously.  Either \code{BIC} (default), or \code{CV}.}
    
  \item{ncv}{The number of cross-validation folds. Default is \code{10}. If \code{method} is not \code{CV}, \code{ncv} is useless.}  
  
  \item{penalty}{The penalty to be applied to the model. Either \code{LASSO} (the default), \code{MCP} or \code{SCAD}.}
  
  \item{lambda}{A user-specified sequence of lambda values.  By default, a sequence of values of length \code{nlam} is computed, equally spaced on the log scale.}
  
  \item{D0}{A user-specified list of initialized values, including five initialized matrix \eqn{S_{(4)}} (called \code{S}), \code{A}, \code{B}, \code{C} and \code{D}. By default, a list of initialization satisfying fixed ranks is computed by random.}
 
  \item{intercept}{A logical value indicating whether the intercept is fitted. Default is \code{TRUE} or set to zero by \code{FALSE}.}
  
  \item{degr}{The number of knots of B-spline base function. Default is \code{3}.}
  
  \item{nlam}{The number of lambda values. Default is \code{20}.}
  
  \item{lam_min}{The smallest value for lambda, as a fraction of lambda.max.  Default is \code{0.01}.}
  
  \item{eps}{Convergence threshhold. The algorithm iterates until the relative change in any coefficient is less than \code{eps}. Default is \code{1e-4}.}
  
  \item{max_step}{Maximum number of iterations. Default is \code{20}.}
              
  \item{eps1}{Convergence threshhold. The Coordinate descent method algorithm iterates until the relative change in any coefficient is less than \code{eps1}. Default is \code{1e-4}.}
  
  \item{max_step1}{The maximum iterates number of coordinate descent method. Default is \code{20}.}
  
  \item{gamma}{The tuning parameter of the MCP or SCAD penalty (see details).}
  
  \item{dfmax}{Upper bound for the number of nonzero coefficients. Default is no upper bound.  However, for large data sets, computational burden may be heavy for models with a large number of nonzero coefficients.}
               
  \item{alpha}{Tuning parameter for the Mnet estimator which controls the relative contributions from the LASSO, MCP or SCAD penalty and the ridge, or L2 penalty.  \code{alpha=1} is equivalent to LASSO, MCP or SCAD penalty, while \code{alpha=0} would be equivalent to ridge regression. However, \code{alpha=0} is not supported; \code{alpha} may be arbitrarily small, but not exactly 0.}
  
}

\details{
  This function gives \code{pq} functional coefficients' estimators of structural marm. A fourth-order tensor with unknown ranks (\eqn{r_1, r_2, r_3, r_4}) need to be estimated. The BIC or CV can be applied to select the optimal regularization parameter, multiple ranks and the number of B-spline basis functions simultaneously. An alternative updating algorithm can be used to update its core tensor and factor matrices based on Tucker decomposition. Group LASSO, SCAD or MCP penalty is applied on the row of the factor matrix \eqn{A} to achieve variable selection.
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
  
  \item{rk_opt}{The optimal ranks and the number of B-spline basis functions that selected by \code{BIC}, or \code{CV}. It is a vector with length 5, which are selected \eqn{r_1}, \eqn{r_2}, \eqn{r_3}, \eqn{r_4} and \eqn{K}.}
  
  \item{lambda.seq}{The sequence of regularization parameter values in the path.}
  
  \item{lambda_opt}{The value of \code{lambda} with the minimum \code{BIC} or \code{CV} value.}
  
  \item{rss}{Residual sum of squares (RSS).}
  
  \item{df}{Degrees of freedom.}
  
  \item{activeX}{The active set of \eqn{X}. A length-\eqn{p} vector.}
  
  \item{opts}{Other related parameters used in algorithm. Some of them are set by default.}
  
  \item{opts}{Other related parameters used in algorithm (especially parameters in peanlty). Some of them are set by default.}

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
  marm4
}

\examples{
library(tensorMARM)
n <- 200; q <- 5; p <- 100; s <- 3; ng = 4
group <- rep(1:ng,each=p/ng)
mydata <- marm4.sim.fbs(n,q,p,s,group)
fit <- with(mydata, marm4.dr(Y,X,group,K,r1_index,r2_index,r3_index,r4_index,D0=D0,nlam=5))
}