\name{fevcov}
\alias{fevcov}
\title{Compute bias corrrected covariance matrix between fixed effects}
\concept{Limited Mobility Bias}
\description{
With a model like 'y = X beta + D theta + F psi + epsilon', where 'D'
and 'F' are matrices with dummy encoded factors, one application of \pkg{lfe} is to
study the variances 'var(D theta), var(F psi)' and covariances 'cov(D theta,
F psi)'
However, if we use
estimates for theta and psi, the resulting variances are biased.
The function \code{fevcov} computes a bias corrected covariance matrix
as described in \cite{Gaure (2014)}.
}

\usage{fevcov(est, alpha=getfe(est), tol=0.01)}

\arguments{
  \item{est}{an object of class '"felm"', the result of a call to
  \code{\link{felm}(keepX=TRUE)}.}

  \item{alpha}{a data frame, the result of a call to \code{\link{getfe}}.}

  \item{tol}{numeric. The absolute tolerance for the bias-corrected correlation.} 
}

\value{
  \code{fevcov} returns a square matrix with the bias corrected covariances. An
  attribute \code{'bias'} contains the biases.
  The bias corrections have been subtracted from the bias estimates.
  I.e. vc = vc' - b, where vc' is the biased variance and b is the bias.
}

\details{
  The \code{tol} argument specifies the tolerance. The tolerance is
  relative for the variances, i.e. the diagonal of the output.
  For the covariances, the tolerance is relative to the square root of
  the product of the variances, i.e. an absolute tolerance for the correlation.
  If a numeric of length  1 \code{tol} specifies the same tolerance for
  all variances/covariances.
  If it is of length 2, \code{tol[1]} specifies the variance tolerance, and
  \code{tol[2]} the covariance tolerance.  \code{tol} can also be a
  square matrix of size \code{length(est$fe)}, in which case the
  tolerance for each variance and covariance is specified individually.

  The function performs no checks for estimability. If the fixed effects
  are not estimable, the result of a call to \code{fevcov} is not
  useable. Moreover, there should be just a single connected component
  among the fixed effects.
}

\seealso{
\code{\link{varvars}}
\code{\link{bccorr}}
}

\note{
Note that if \code{est} is the result of a call to \code{\link{felm}}
with \code{keepX=FALSE} (the default), the biases will be computed
as if the covariates X are independent of the factors. This will be
faster (typically by a factor of approx. 4), and possibly wronger.
Note also that the computations performed by this function are
non-trivial, they may take quite some time.  It would be wise to start
out with quite liberal tolerances, e.g. \cite{tol=0.1}, to
get an idea of the time requirements.
Currently, only i.i.d. residuals are supported for this function.
}


\examples{
x <- rnorm(500)
x2 <- rnorm(length(x))

## create individual and firm
id <- factor(sample(40,length(x),replace=TRUE))
firm <- factor(sample(30,length(x),replace=TRUE,prob=c(2,rep(1,29))))
foo <- factor(sample(20,length(x),replace=TRUE))
## effects
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))
foo.eff <- rnorm(nlevels(foo))
## left hand side
id.m <- id.eff[id]
firm.m <- 2*firm.eff[firm]
foo.m <- 3*foo.eff[foo]
y <- x + 0.25*x2 + id.m + firm.m + foo.m + rnorm(length(x))

# make a data frame
fr <- data.frame(y,x,x2,id,firm,foo)
## estimate and print result
est <- felm(y ~ x+x2|id+firm+foo, data=fr, keepX=TRUE)
# find bias corrections
fevcov(est)
}

\references{
  Gaure, S. (2014), \cite{Correlation bias correction in two-way
  fixed-effects linear regression}, Stat 3(1):379-390, 2014.
}
