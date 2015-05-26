\name{condfstat}
\alias{condfstat}
\title{Compute conditional F statistic for weak instruments in an
  IV-estimation with multiple endogenous variables}
\description{
  When using multiple instruments for multiple endogenous variables,
  the ordinary individual t-tests for the instruments in the first stage
  do not always reveal a weak set of instruments.
  Conditional F statistics can be used for such testing.
}
\usage{
condfstat(object, type='default', quantiles=0.0, bN=100L)
}

\arguments{
  \item{object}{object of class \code{"felm"}, a result of a call
  to \code{\link{felm}}.}
  \item{type}{character. Error structure. Passed to
  \code{\link{waldtest}}. If \code{NULL}, both iid and robust Fs are returned.}
  \item{quantiles}{numeric. Quantiles for bootstrap.}
  \item{bN}{integer. Number of bootstrap samples.}
}

\value{
A p x k matrix, where k is the number of endogenous
variables. Each row are the conditional F statistics
on a residual equation as described
in \cite{Sanderson and Windmeijer (2014)}, for a certain
error structure.  The default is to use iid, or cluster if a cluster
was specified to \code{\link{felm}}. The third choice is
\code{'robust'}, for heteroskedastic errors.
If \code{type=NULL}, iid and robust Fs are returned, and cluster, if
that was specified to \code{felm}.

Note that for these F statistics it is not the p-value that matters, it
is the F statistic itself which (coincidentally) pops up in the denominator for the
asymptotic bias of the IV estimates, and thus a large F is beneficial.

}

\details{

  IV coefficient estimates are not normally distributed, in particular
  they do not have the right expectation.  They follow a quite
  complicated distribution which is fairly close to normal if the
  instruments are good.  The conditional F-statistic is a measure of how
  good the instruments are.  If the F is large, the instruments are
  good, and any bias due to the instruments is small compared to the
  estimated standard errors, and also small relative to the bias in
  OLS. See \cite{Sanderson and Windmeijer (2014)} and \cite{Stock and
    Yogo (2004)}.  If F is small, the bias can be large compared to the
  standard error.


  If \code{any(quantiles > 0.0)}, a bootstrap with \code{bN} samples
  will be performed to estimate quantiles of the endogenous parameters which
  includes the variance both from the 1st and 2nd stage.  The result is
  returned in an array attribute \code{quantiles} of the value returned
  by \code{condfstat}.  The argument \code{quantiles} can be a vector
  to estimate more than one quantile at once.  If \code{quantiles=NULL},
  the bootstrapped estimates themselves are returned.  The bootstrap is
  normally much faster than running \code{felm} over and over
  again. This is so because all exogenous variables are projected out of
  the equations before doing the bootstrap.

}

\note{
Please note that \code{condfstat} does not work with the old syntax for IV
in \code{\link{felm}(...,iv=)}. The new multipart syntax must be used.
}

\references{
Sanderson, E. and F. Windmeijer (2014)
\cite{A weak instrument F-test in linear IV models with multiple
endogenous variables}, Disc. Paper 14/644, Univ of Bristol.
\url{http://www.efm.bris.ac.uk/economics/working_papers/pdffiles/dp14644.pdf}

Stock, J.H. and M. Yogo (2004)
\cite{Testing for weak instruments in linear IV regression},
\url{http://ssrn.com/abstract=1734933} in
\cite{Identification and inference for econometric models: Essays in honor of
Thomas Rothenberg}, 2005.
}

\examples{
z1 <- rnorm(4000)
z2 <- rnorm(length(z1))
u <- rnorm(length(z1))
# make x1, x2 correlated with errors u

x1 <- z1 + z2 + 0.2*u + rnorm(length(z1))
x2 <- z1 + 0.94*z2 - 0.3*u + rnorm(length(z1))
y <- x1 + x2 + u
est <- felm(y ~ 1 | 0 | (x1 | x2 ~ z1 + z2))
summary(est)
\dontrun{
summary(est$stage1, lhs='x1')
summary(est$stage1, lhs='x2')
}

# the joint significance of the instruments in both the first stages are ok:
t(sapply(est$stage1$lhs, function(lh) waldtest(est$stage1, ~z1|z2, lhs=lh)))
# everything above looks fine, t-tests for instruments, 
# as well as F-tests for excluded instruments in the 1st stages.
# The conditional F-test reveals that the instruments are jointly weak
# (it's close to being only one instrument, z1+z2, for both x1 and x2)
condfstat(est, quantiles=c(0.05, 0.95))
}

