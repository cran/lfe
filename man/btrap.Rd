\name{btrap}
\alias{btrap}
\title{Bootstrap standard errors for the group fixed effects}
\concept{Bootstrap}
\description{
  Bootstrap standard errors for the group fixed effects which were swept out during
  an estimation with \code{\link{felm}}.
}
\usage{
btrap(alpha,obj,N=100,ef=NULL,eps=getOption('lfe.eps'),threads=getOption('lfe.threads'))
}


\arguments{
\item{alpha}{data frame returned from \code{\link{getfe}}}
\item{obj}{object of class \code{"felm"}, usually, a result of a call to
  \code{\link{felm}}}
\item{N}{integer.  The number of bootstrap iterations}
\item{ef}{function.  An estimable function such as in
   \code{\link{getfe}}. The default is to use the one used on \code{alpha}}
\item{eps}{double. Tolerance for centering, as in getfe}
\item{threads}{integer.  The number of threads to use}
}

\value{
A data-frame of the same size as alpha is returned, with standard errors
filled in.
}

\details{
  The bootstrapping is done in parallel if \code{threads > 1}.  However,
  for better utilization of parallelism one should load one of the
  packages \pkg{doMC} or \pkg{doSMP}.  The former seems to perform
  better in this context.
  \code{\link{btrap}} is run automatically from \code{\link{getfe}} if
  \code{se=TRUE} is specified.
}

\examples{
## create covariates
x <- rnorm(5000)
x2 <- rnorm(length(x))

## create individual and firm
id <- factor(sample(1000,length(x),replace=TRUE))
firm <- factor(sample(500,length(x),replace=TRUE))

## effects
id.eff <- rlnorm(nlevels(id))
firm.eff <- rexp(nlevels(firm))

## left hand side
y <- x + 0.25*x2 + id.eff[id] + firm.eff[firm] + rnorm(length(x))

## estimate and print result
est <- felm(y ~ x+x2 + G(id)+G(firm))
summary(est)
## extract the group effects
alpha <- getfe(est)
head(alpha)
## bootstrap standard errors
head(btrap(alpha,est))

## bootstrap som differences
ef <- function(v) {
  w <- c(v[2]-v[1],v[3]-v[2],v[3]-v[1])
  names(w) <-c('id2-id1','id3-id2','id3-id1')
  w
}
attr(ef,'noP') <- TRUE
head(btrap(alpha,est,ef=ef))

}
