\name{getfe}
\alias{getfe}
\title{Retrieve the group fixed effects}
\concept{Kaczmarz Method}
\description{
  Compute the group fixed effects which were swept out during
  an estimation with \code{\link{felm}}.
}
\usage{
getfe(obj,references=NULL,se=FALSE,method='kaczmarz',ef=defaultef(obj),bN=100) 
}

\arguments{
\item{obj}{object of class \code{"felm"}, usually, a result of a call to
  \code{\link{felm}}}
\item{references}{a vector of strings.  If there are more than two
  factors and you have prior knowledge of what the reference levels
  should be, use this argument, like \code{references='id.23'}.  Not
  used with \code{method='kaczmarz'}}
\item{se}{logical.  Set to TRUE if standard errors for the group
  effects are wanted.  This is \strong{very} time-consuming for large problems, so
  leave it as FALSE unless absolutely needed.}
\item{method}{character string.  Either 'cholesky', the default, or
  'kaczmarz'.  The latter is often very fast and consumes little memory, it will yield the smallest norm
  solution.
  }
\item{ef}{function. A function of one variable, a vector of group fixed
  effects. This function should be estimable and is used to transform
  the raw-coefficients from the kaczmarz-method.  The default value is
  one which picks one reference in each component.  Can be set to NULL
  to yield the minimal-norm solution from the kaczmarz-method.  Other
  functions may be added later, upon request.}
\item{bN}{integer.  The number of bootstrap runs when standard errors
  are requested}
}

\value{
The function \code{getfe} computes and returns a data frame containing
the group fixed effects.  It has the columns
\code{c('effect','se','obs','comp','fe','idx')}

\itemize{
  \item \code{effect} is the estimated effect.
  \item \code{se} is the standard error.
  \item \code{obs} is the number of observations of this level.
  \item \code{comp} is the graph-theoretic component number, useful for
  interpreting the effects.
  \item \code{fe} is the name of factor.
  \item \code{idx} is the level of the factor.
}

With the Kaczmarz-method it's possible to specify a different estimable
function. If the function returns a vector of a different length, or the
function has the attribute 'noP' set to TRUE, only the first two columns
are returned.
}

\details{
  For the case with two factors (the \code{fl}-argument to
  \code{\link{felm}}), one reference in each connection component is
  adequate when interpreting the results.

  For three or more factors, no such easy method is known;
  for the \code{"cholesky"} method- reference levels are found by analyzing the pivoted
  Cholesky-decomposition of a slightly perturbed system.  The
  \code{"kaczmarz"} method provides no rank-deficiency analysis, it is
  assumed that the factors beyond the two first contribute nothing to
  the rank-deficiency, so one reference in each is used.

  If there are more than two factors, only the first two will be used
  to report connection components.  In this case, it is not known
  which graph theoretic concept may be used to analyze the rank-deficiency.

  The standard errors returned by the Kaczmarz-method are bootstrapped,
  keeping the other coefficients (from \code{\link{felm}}) constant,
  i.e. they are from the variance when resampling the
  residuals. 
}

\examples{
## create covariates
x <- rnorm(50000)
x2 <- rnorm(length(x))

## create individual and firm
id <- factor(sample(5000,length(x),replace=TRUE))
firm <- factor(sample(3000,length(x),replace=TRUE))

## effects
id.eff <- rlnorm(nlevels(id))
firm.eff <- rexp(nlevels(firm))

## left hand side
y <- x + 0.25*x2 + id.eff[id] + firm.eff[firm] + rnorm(length(x))

## estimate and print result
est <- felm(y ~ x+x2 + G(id)+G(firm))
summary(est)
## extract the group effects
alpha <- getfe(est,se=TRUE)

## find some differences, with standard errors
head(alpha,3)
ef <- function(v) {
  w <- c(v[2]-v[1],v[3]-v[2],v[3]-v[1])
  names(w) <-c('id2-id1','id3-id2','id3-id1')
  w
}
attr(ef,'noP') <- TRUE
getfe(est,ef=ef,se=TRUE)

\dontrun{
summary(lm(y ~ x+x2+id+firm-1))
}
}
\keyword{regression}
\keyword{models}
