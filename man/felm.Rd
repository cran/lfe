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

\usage{ felm(model,fl,data) }

\arguments{
  \item{model}{an object of class '"formula"' (or one that can be
  coerced to that class: a symbolic description of the model to be
  fitted. Similarly to 'lm'}
  \item{fl}{a list of factors defining the groups. The list
  must have unique names, see Examples section.}
  \item{data}{a data frame containing the variables of the model}
}

\value{
  \code{felm} return an object of \code{class} \code{"lm.fe"}.  The
  content of this is currently not documented, but there is a method
  \code{\link{print.lm.fe}} which may be used to display the content in
  a similar way as for \code{\link{lm}}.
}

\seealso{\code{\link{getfe}}}
\examples{
## create covariates
x <- rnorm(100)
x2 <- rnorm(length(x))

## individual and firm
id <- factor(sample(10,length(x),replace=TRUE))
firm <- factor(sample(3,length(x),replace=TRUE))

## effects for them
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

## left hand side
y <- x + 0.5*x2 + id.eff[id] + firm.eff[firm] + rnorm(length(x))

## estimate and print result
est <- felm(y ~ x+x2,fl=list(id=id,firm=firm),data=data.frame(y,x,x2))
print(est)
\dontrun{print(getfe(est))}
}
\keyword{regression}
\keyword{models}
