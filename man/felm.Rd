\name{felm}
\alias{felm}
\title{Fitting linear models with multiple group fixed effects}
\description{
'felm' is used to fit linear models with multiple group fixed effects,
similarly to lm.  It uses the Method of Alternating projections to
sweep out multiple group effects from the normal equations.

This function is intended for use with large datasets with multiple group
effects of large cardinality.  If dummy-encoding the group effects
results in a manageable number of coefficients, you are probably better
off by using \code{\link{lm}} or package \pkg{lme4}.
}

\usage{ felm(formula,fl,data) }

\arguments{
  \item{formula}{an object of class '"formula"' (or one that can be
  coerced to that class: a symbolic description of the model to be
  fitted. Similarly to 'lm'.  Grouping factors are coded as G(f) (with a
  capital \code{G}).}

\item{fl}{a list of factors defining the groups. An alternative to the
   G(f) syntax for backward compatibility.  The list must have unique
   names, see Examples section.}

  \item{data}{a data frame containing the variables of the model}
}

\value{
  \code{felm} return an object of \code{class} \code{"felm"}.  The
  content of this is currently not documented, but the generic
  \code{summary}-method will yield a summary which may be \code{print}'ed.
}

\note{
The standard errors are adjusted for the reduced degrees of freedom
coming from the dummies which are implicitly present.  In the case of
two factors, the exact number of dummies is fast to compute.  If there
are more factors, the number of dummies may be slightly over-estimated,
leading to slightly too large standard errors.
}

\seealso{\code{\link{getfe}}}
\examples{
## create covariates
x <- rnorm(100)
x2 <- rnorm(length(x))

## individual and firm
id <- factor(sample(10,length(x),replace=TRUE))
firm <- factor(sample(3,length(x),replace=TRUE,prob=c(2,1,1)))

## effects for them
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

## left hand side
y <- x + 0.5*x2 + id.eff[id] + firm.eff[firm] + rnorm(length(x))

## estimate and print result
est <- felm(y ~ x+x2+G(id)+G(firm))
summary(est)
## compare with lm
summary(lm(y ~ x + x2 + id + firm-1))
## alternatively
\dontrun{felm(y ~ x + x2,fl=list(id=id,firm=firm))}

\dontrun{print(getfe(est))}

}
\keyword{regression}
\keyword{models}
