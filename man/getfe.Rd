\name{getfe}
\alias{getfe}
\title{Retrieve the group fixed effects}
\description{
  Compute the group fixed effects which were swept out during
  an estimation with \code{\link{felm}}.
}
\usage{
getfe(obj,data) 
}

\arguments{
\item{obj}{object of class \code{"fe.lm"}, usually, a result of a call to
  \code{\link{felm}}}
\item{data}{a data frame}
}

\value{
The function \code{getfe} computes and returns a data frame containing
the group fixed effects.
It has the columns \code{c('effect','obs','comp','fe','idx')}

\itemize{
  \item \code{effect} is the estimated effect.
  \item \code{obs} is the number of observations of this level.
  \item \code{comp} is the graph-theoretic component number, useful for
  interpreting the effects.
  \item \code{fe} is the name of factor.
  \item \code{idx} is the level of the factor.
}
}

\examples{
## create covariates
x <- rnorm(100000)
x2 <- rnorm(length(x))

## create individual and firm
id <- factor(sample(20,length(x),replace=TRUE))
firm <- factor(sample(8,length(x),replace=TRUE,prob=c(2,rep(1,7))))

## effects
id.eff <- rlnorm(nlevels(id))
firm.eff <- rexp(nlevels(firm))

## left hand side
y <- x + 0.25*x2 + id.eff[id] + firm.eff[firm] + rnorm(length(x))

## make a data frame
dset <- data.frame(y,x,x2,id,firm)

## estimate and print result
est <- felm(y ~ x+x2,fl=list(id=id,firm=firm),data=dset)
print(est)

## extract the group effects, and print them
print(getfe(est))

## compare with an ordinary lm
summary(lm(y ~ x+x2+id+firm-1,data=dset))

## compare with lmer
if(require(lme4))
coef(lmer(y ~ x + x2 + (1|id) + (1|firm),data=dset))
}
\keyword{regression}
\keyword{models}
