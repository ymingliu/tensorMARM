\name{tensorMARM-package}
\alias{tensorMARM-package}
\alias{tensorMARM}
\docType{package}
\title{
  Multivariate additive regression for multi-view data via tensor estimation
}
\description{
  For an integrative multi-view multivariate additive model (MARM), the B-splines are applied to approximate the component functions. We treat the coefficients as multiple third-order tensors or even a fourth-order tensor in some specific conditions (called structural MARM). With the tensor low-rankness, the Tucker decomposition and group sparse penalty (lasso, mcp or scad) reduce the number of parameters. An alternative updating algorithm based on the coordinate descent strategy is used to estimate the core tensors and factor matrices, and further additive functions.

}
\details{
  This pacakge includes four main functions and four generating functions. \code{marm3} and \code{marm3.dr} yield the estimator of MARM. The difference is the former requires the fixed ranks and the latter can search the optimal ranks and regularization parameter simultaneously by BIC or CV method. \code{marm3.sim.fbs} and \code{marm3.sim.fsin} generate data of scenario I and II respectively. Scenario I assumes that the true functions are exactly residing in the space of B-spline basis functions. Scenario II assumes that the true functions are some linear combination of \eqn{sin(2\pi x)} and \eqn{cos(\pi x)}. Similarly, \code{marm4} and \code{marm4.dr} yield the estimator of structural MARM. \code{marm4.sim.fbs} and \code{marm4.sim.fsin} are two generating functions of scenario I and II, which have the same assumptions as MARM.
}
\author{
Xu Liu

Maintainer: Xu Liu <liu.xu@sufe.edu.cn>
}
\references{
  Multivariate additive regression for multi-view data via tensor estimation.
}
\keyword{ Group sparsity; Tensor low-rankness; Tucker decomposition; Multivariate additive regression; Multi-view data. }
