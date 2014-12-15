\name{bccorr}
\alias{bccorr}
\title{Compute bias corrrected correlation between fixed effects}
\concept{Limited Mobility Bias}
\description{
With a model like 'y = X beta + D theta + F psi + epsilon', where 'D'
and 'F' are matrices with dummy encoded factors, one application of \pkg{lfe} is to
study the correlation 'cor(D theta, F psi)'.  However, if we use
estimates for theta and psi, the resulting correlation is biased.
The function \code{bccorr} computes a bias corrected correlation
as described in \cite{Gaure (2014)}.
}

\usage{bccorr(est, alpha=getfe(est), corrfactors=1L:2L,
                  nocovar=is.null(est$X) && length(est$fe)==2,
                  tol=0.01, maxsamples=Inf)}

\arguments{
  \item{est}{an object of class '"felm"', the result of a call to
  \code{\link{felm}(keepX=TRUE)}.}

  \item{alpha}{a data frame, the result of a call to \code{\link{getfe}}.}

  \item{corrfactors}{integer or character vector of length 2. The factors to
    correlate. The default is fine if there are only two factors in the model.}
  
  \item{nocovar}{logical. Assume no other covariates than the two
  factors are present, or that they are uncorrelated with them.}

  \item{tol}{The absolute tolerance for the bias-corrected correlation.} 

  \item{maxsamples}{Maximum number of samples for the trace sample means estimates}
}

\value{
  \code{bccorr} returns a named integer vector with the following fields:

  \item{corr}{the bias corrected correlation.}
  \item{v1}{the bias corrected variance for the first factor specified
  by \code{corrfactors}.}
  \item{v2}{the bias corrected variance for the second factor.}
  \item{cov}{the bias corrected covariance between the two factors.}
  \item{d1}{the bias correction for the first factor.}
  \item{d2}{the bias correction for the second factor.}
  \item{d12}{the bias correction for covariance.}

  The bias corrections have been subtracted from the bias estimates.
  E.g. v2 = v2' - d2, where v2' is the biased variance.
}

\details{
The bias expressions from \cite{Andrews et al.} are of the form 'tr(A B^\{-1\} C)'
where A, B, and C are  matrices too large to be handled
directly. \code{bccorr} estimates the trace by using the formula \eqn{tr(M) = E(x^t M x)}
where x is a vector with coordinates drawn uniformly from the set \{-1,1\}.
More specifically, the expectation is estimated by
sample means, i.e. in each sample a vector x is drawn, the
equation 'Bv = Cx' is solved by a conjugate gradient method, and the
real number \eqn{x^t Av} is computed. 

There are three bias corrections, for the variances of D theta (vD) and
F psi (vF), and their covariance (vDF).The correlation is computed as
\code{rho <- vDF/sqrt(vD*vF)}.  The variances are estimated to a
relative tolerance specified by the argument \code{tol}. The covariance
bias is estimated to an absolute tolerance in the correlation \code{rho}
(conditional on the already bias corrected vD and vF) specified by
\code{tol}.  The CG algortithm does not need to be exceedingly precise,
it is terminated when the solution reaches a precision which is
sufficient for the chosen precision in \code{vD, vF, vDF}.
}

\note{
Note that if \code{est} is the result of a call to \code{\link{felm}}
with \code{keepX=FALSE} (the default), the correlation will be computed
as if the covariates X are independent of the two factors. This will be
faster (typically by a factor of approx. 4), and possibly wronger.

Note also that the computations performed by this function are
non-trivial, they may take quite some time.  It would be wise to start
out with quite liberal tolerances, e.g. \cite{tol=0.1}, to
get an idea of the time requirements.

The algorithm used is not very well suited for small datasets with only
a few thousand levels in the factors.
}

\references{
  Gaure, S. (2014), \cite{Correlation bias correction in two-way
  fixed-effects linear regression}, Stat 3(1), 2014, 379-390.
}
