% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bccorr.R
\name{varvars}
\alias{varvars}
\title{Compute the variance of the fixed effect variance estimate}
\usage{
varvars(est, alpha = getfe(est), tol = 0.01, biascorrect = FALSE, lhs = NULL)
}
\arguments{
\item{est}{an object of class '"felm"', the result of a call to
\verb{[felm](keepX=TRUE)}.}

\item{alpha}{a data frame, the result of a call to \code{\link[=getfe]{getfe()}}.}

\item{tol}{numeric. The absolute tolerance for the bias-corrected
correlation.}

\item{biascorrect}{logical. Should the estimates be bias corrected?}

\item{lhs}{character. Name of left hand side if multiple left hand sides.}
}
\value{
\code{varvars} returns a vector with a variance estimate for each
fixed effect variance.  I.e. for the diagonal returned by
\code{\link[=fevcov]{fevcov()}}.
}
\description{
Compute the variance of the fixed effect variance estimate
}
\details{
With a model like \eqn{y = X\beta + D\theta + F\psi + \epsilon}, where \eqn{D} and
\eqn{F} are matrices with dummy encoded factors, one application of \pkg{lfe} is
to study the variances \eqn{var(D\theta)}, \eqn{var(F\psi)} and covariances
\eqn{cov(D\theta, F\psi)}. The function \code{\link[=fevcov]{fevcov()}} computes bias corrected
variances and covariances.  However, these variance estimates are still
random variables for which \code{\link[=fevcov]{fevcov()}} only estimate the
expectation. The function \code{varvars} estimates the variance of these
estimates.

This function returns valid results only for normally distributed residuals.
Note that the estimates for the fixed effect variances from
\code{\link[=fevcov]{fevcov()}} are not normally distributed, but a sum of chi-square
distributions which depends on the eigenvalues of certain large matrices. We
do not compute that distribution. The variances returned by \code{varvars}
can therefore \emph{not} be used directly to estimate confidence intervals,
other than through coarse methods like the Chebyshev inequality. These
estimates only serve as a rough guideline as to how wrong the variance
estimates from \code{\link[=fevcov]{fevcov()}} might be.

Like the fixed effect variances themselves, their variances are also biased
upwards.  Correcting this bias can be costly, and is therefore by default
switched off.

The variances tend to zero with increasing number of observations. Thus, for
large datasets they will be quite small.
}
\note{
The \code{tol} argument specifies the tolerance as in
\code{\link[=fevcov]{fevcov()}}.  Note that if \code{est} is the result of a call to
\code{\link[=felm]{felm()}} with \code{keepX=FALSE} (the default), the variances will
be estimated as if the covariates X are independent of the factors.  There
is currently no function available for estimating the variance of the
covariance estimates from \code{\link[=fevcov]{fevcov()}}.

The cited paper does not contain the expressions for the variances computed
by \code{varvars} (there's a 10 page limit in that journal), though they can
be derived in the same fashion as in the paper, with the formula for the
variance of a quadratic form.
}
\examples{

x <- rnorm(500)
x2 <- rnorm(length(x))

## create individual and firm
id <- factor(sample(40, length(x), replace = TRUE))
firm <- factor(sample(30, length(x), replace = TRUE, prob = c(2, rep(1, 29))))
foo <- factor(sample(20, length(x), replace = TRUE))
## effects
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))
foo.eff <- rnorm(nlevels(foo))
## left hand side
id.m <- id.eff[id]
firm.m <- 2 * firm.eff[firm]
foo.m <- 3 * foo.eff[foo]
y <- x + 0.25 * x2 + id.m + firm.m + foo.m + rnorm(length(x))

# make a data frame
fr <- data.frame(y, x, x2, id, firm, foo)
## estimate and print result
est <- felm(y ~ x + x2 | id + firm + foo, data = fr, keepX = TRUE)
alpha <- getfe(est)
# estimate the covariance matrix of the fixed effects
fevcov(est, alpha)
# estimate variances of the diagonal
varvars(est, alpha)

}
\references{
Gaure, S. (2014), \cite{Correlation bias correction in two-way
fixed-effects linear regression}, Stat 3(1):379-390, 2014.
}
\seealso{
\code{\link[=bccorr]{bccorr()}} \code{\link[=fevcov]{fevcov()}}
}
