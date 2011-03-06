\name{getfe}
\alias{getfe}
\title{Retrieve the group fixed effects}
\description{
  Compute the group fixed effects which were swept out during
  an estimation with \code{\link{felm}}.
}
\usage{
getfe(obj,references=NULL,se=FALSE) 
}

\arguments{
\item{obj}{object of class \code{"felm"}, usually, a result of a call to
  \code{\link{felm}}}
\item{references}{a vector of strings.  If there are more than two
  factors and you have prior knowledge of what the reference levels
  should be, use this argument, like \code{references='id.23'}.}
\item{se}{logical.  Set to TRUE if standard errors for the group
  effects are wanted.  This is *very* time-consuming for large problems, so
  leave it as FALSE unless absolutely needed.}
}

\value{
The function \code{getfe} computes and returns a data frame containing
the group fixed effects.  It has the columns
\code{c('effect','obs','comp','fe','idx')}

\itemize{
  \item \code{effect} is the estimated effect.
  \item \code{obs} is the number of observations of this level.
  \item \code{comp} is the graph-theoretic component number, useful for
  interpreting the effects.
  \item \code{fe} is the name of factor.
  \item \code{idx} is the level of the factor.
}
}

\details{
  For the case with two factors (the \code{fl}-argument to
  \code{\link{felm}}), one reference in each connection component is
  chosen automatically.  For three or more factors, no such easy method is known,
  reference levels are found by analyzing the pivoted
  Cholesky-decomposition of a slightly perturbed system.
}

\examples{
## create covariates
x <- rnorm(100000)
x2 <- rnorm(length(x))

## create individual and firm
id <- factor(sample(3000,length(x),replace=TRUE))
firm <- factor(sample(1500,length(x),replace=TRUE,prob=c(2,rep(1,1499))))

## effects
id.eff <- rlnorm(nlevels(id))
firm.eff <- rexp(nlevels(firm))

## left hand side
y <- x + 0.25*x2 + id.eff[id] + firm.eff[firm] + rnorm(length(x))

## estimate and print result
cat(date(),'Start lfe estimation\n')
est <- felm(y ~ x+x2 + G(id)+G(firm))
summary(est)
## extract the group effects
alpha <- getfe(est)
cat(date(),'Done\n')
\dontrun{
cat(date(),'Compare with lm estimation\n')
lmsum <- summary(lm(y ~ x+x2+id+firm-1))
cat(date(),'Done\n')
}
}
\keyword{regression}
\keyword{models}
