\name{SurvGPR}
\alias{SurvGPR}
\title{Fit a Gaussian process regression model to right-censored survival time data.}
\usage{
SurvGPR(time, status, Z, K, tol = 1e-7, max.iter = 100, max.iter.MM = 100, quiet = FALSE, 
        max.samples = 1e5, kern.type = c("K+I", "multi-K"))
}
\arguments{
 \item{time}{An \eqn{n}-variate or containing the failure/censoring times (on the original scale -- NOT log-transformed).}
  \item{status}{An \eqn{n}-variate binary vector of same length as time -- 0 indicates censored, 1 indicates failure. }
  \item{Z}{An \eqn{n \times q} design matrix for the linear mean function. Note that the first column should contain all ones. We recommend construction using \code{model.matrix}. }
  \item{K}{Candidate kernel matrices in the form of an array of dimension \eqn{n \times n \times M}. }
  \item{tol}{The convergence tolerance. Default is \code{1e-7}. }
  \item{max.iter.MM}{The maximum number of iterations for the inner M-step algorithm. }
  \item{max.iter}{The maximum number of total EM-iterations. }
  \item{kern.type}{A character argument -- either \code{K+I} when \eqn{M=1} or \code{multi-K} when \eqn{M>1}. Incorrect compatability with input \eqn{K} will produce an error. }
  \item{quiet}{\code{TRUE/FALSE} -- print algorithm progress?}
  \item{max.samples}{An upper bound on \eqn{s_k}, the Monte-Carlo sample size for the \eqn{k}th iteration. Note that the final imputed values of log-survival for censored subjects will be the average of \code{max.samples} Monte-Carlo draws. }
}
\value{
   \item{beta}{\eqn{\hat{\beta}}: The estimated regression coefficient vector corresponding to the columns of \code{Z}. }
  \item{sigma2}{\eqn{\hat{\sigma}^2}: The estimated variance components -- a vector of length \eqn{M+1}, with the final element corresponding to the variance of \eqn{\epsilon}. }
  \item{Tout}{The log-failure and imputed log-failure times obtained from our MCEM algorithm. These are primarily to be used in the prediction function. }
}

\description{
  A function for fitting a Gaussian process regression model to right-censored survival time data. 
}
