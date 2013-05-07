\name{felm}
\alias{felm}
\title{Fitting linear models with multiple group fixed effects}
\description{
'felm' is used to fit linear models with multiple group fixed effects,
similarly to lm.  It uses the Method of Alternating projections to
sweep out multiple group effects from the normal equations before
estimating the remaining coefficients with OLS.

This function is intended for use with large datasets with multiple group
effects of large cardinality.  If dummy-encoding the group effects
results in a manageable number of coefficients, you are probably better
off by using \code{\link{lm}}.
}

\usage{ felm(formula, data, iv=NULL, clustervar=NULL, exactDOF=FALSE) }

\arguments{
  \item{formula}{an object of class '"formula"' (or one that can be
  coerced to that class: a symbolic description of the model to be
  fitted. Similarly to 'lm'.  Grouping factors are coded as G(f) (with a
  capital \code{G}).}

  \item{data}{a data frame containing the variables of the model}

  \item{iv}{a formula describing an instrumented variable. Estimated via two step OLS}
  
  \item{clustervar}{a string or factor.  Either the name of a variable or a factor. Used for
        computing clustered standard errors.}
  \item{exactDOF}{logical. If more than two factors, the degrees of
        freedom used to scale the covariance matrix (and the standard
        errors) is normally estimated. Setting \code{exactDOF=TRUE}
        causes \code{felm} to compute it, but this may fail if there are
        too many levels in the factors.}
}

\value{
  \code{felm} returns an object of \code{class} \code{"felm"}.  The
  content of this is currently not documented, but the generic
  \code{summary}-method will yield a summary which may be \code{print}'ed.
  The object has some resemblance to the an \code{lm} object, and some
  postprocessing methods designed for \code{lm} may happen to work.
}

\note{
The standard errors are adjusted for the reduced degrees of freedom
coming from the dummies which are implicitly present.  In the case of
two factors, the exact number of implicit dummies is easy to compute.  If there
are more factors, the number of dummies is estimated by assuming there's
one reference-level for each factor, this may be a slight over-estimation,
leading to slightly too large standard errors. Setting \code{exactDOF}
computes the exact degrees of freedom with methods in package \pkg{Matrix}.

For the IV-option, it is only necessary to include the instruments on the
right hand side.  The other covariates, from \code{formula}, are added automatically
in the first step.  See the examples. \code{iv} can also be a list of
formulas if more than one variable is instrumented.

Ideally, the \code{clustervar} should have been an option to the \code{summary}-function instead.
However, this would require keeping a copy of the data-matrix in the returned
structure.  Since this function is intended for very large datasets, we discard
the data-matrix to save memory, keeping only residuals and other summary statistics.

Note that the syntax of the \code{felm}-function has changed, it does no longer
allow a separate specification of the group factors, they must be specified
with the \code{G()}-syntax.  The old felm is still available as \code{lfe:::felm.old}, but
it will no longer be maintained.
}

\seealso{\code{\link{getfe}}}
\examples{
## create covariates
x <- rnorm(1000)
x2 <- rnorm(length(x))

## individual and firm
id <- factor(sample(20,length(x),replace=TRUE))
firm <- factor(sample(13,length(x),replace=TRUE))

## effects for them
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

## left hand side
u <- rnorm(length(x))
y <- x + 0.5*x2 + id.eff[id] + firm.eff[firm] + u

## estimate and print result
est <- felm(y ~ x+x2+G(id)+G(firm))
summary(est)
## compare with lm
summary(lm(y ~ x + x2 + id + firm-1))
## alternatively
\dontrun{felm(y ~ x + x2,fl=list(id=id,firm=firm))
  getfe(est)
}

## make an iv-example, Q is instrumented by x3, report robust s.e.
x3 <- rnorm(length(x))
Q <- 0.3*x3 + x + 0.2*x2 + id.eff[id] + 0.7*u + rnorm(length(x),sd=0.3)
y <- y + Q
ivest <- felm(y ~ x + x2 + G(id)+G(firm) + Q, iv=Q ~ x3)
summary(ivest,robust=TRUE)
}
\keyword{regression}
\keyword{models}
