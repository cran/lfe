\name{condfstat}
\alias{condfstat}
\title{Compute conditional F statistic for weak instruments in an
  IV-estimation with multiple endogenous variables.}
\description{
  When using multiple instruments for multiple endogenous variables,
  the ordinary individual t-tests for the instruments in the first stage
  do not always reveal a weak set of instruments.
  Conditional F statistics can be used for such testing.
}
\usage{
condfstat(object, type='default')
}

\arguments{
  \item{object}{object of class \code{"felm"}, a result of a call
  to \code{\link{felm}}.}
  \item{type}{character. Error structure. Passed to
  \code{\link{waldtest}}. If \code{NULL}, both iid and robust Fs are returned.}
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
This F statistic is not known to have the same interpretation for
heteroskedastic and clustered errors.
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

}

\examples{
set.seed(43)
z1 <- rnorm(4000)
z2 <- rnorm(length(z1))
u <- rnorm(length(z1))
# make x1, x2 correlated with errors u

x1 <- z1 + z2 + 0.2*u + rnorm(length(z1))
x2 <- z1 + 0.95*z2 - 0.3*u + rnorm(length(z1))
y <- x1 + x2 + u
est <- felm(y ~ 1 | 0 | (x1 | x2 ~ z1 + z2))
\dontrun{
summary(est$stage1, lhs='x1')
summary(est$stage1, lhs='x2')
}
summary(est)
# everything above looks fine, t-tests for instruments, 
# as well as F-tests for excluded instruments in the 1st stages.
# The standard errors are large, though
# However, the conditional F-test reveals that the instruments are weak
# (it's close to being only one instrument, z1+z2, for both x1 and x2)
condfstat(est,type=NULL)
}

