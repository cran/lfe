\name{lfe-package}
\alias{lfe-package}
\alias{lfe}
\docType{package}
\title{Linear Group Fixed Effects}
\concept{Method of Alternating Projections}
\description{
The package uses the Method of Alternating Projections to estimate
linear models with multiple group fixed effects.
}

\details{

This package is intended for linear models with multiple group fixed effects.
It performs no other functions than \code{\link[stats]{lm}} or package
\pkg{lme4}, but it uses a special method for projecting out multiple group fixed
effects from the normal equations, hence it is faster. It is a
generalization of the within-groups estimator.  This may be required if
the groups have high cardinality (many levels), resulting in tens
or hundreds of thousands of dummy-variables.  It is also useful if one
only wants to control for the group effects, without actually computing
them.  The package is not able to compute standard errors for the group
effects.

The estimation is done in two steps.  First the other coefficients are
estimated with the function \code{\link{felm}} by centering on all the
group means.  Then the group effects are extracted (if needed)
with the function \code{\link{getfe}}.  There's also a function
\code{\link{demeanlist}} which just does the centering on an arbitrary
matrix, and there's a function \code{\link{compfactor}} which computes
the connection components (which is used for interpreting the group
effects when there are only two factors, see the Abowd et al
references), they are also returned by \code{\link{getfe}}).

The centering on the means is done with a tolerance.  This tolerance is
set by \code{options(lfe.eps=1e-7)}, its default is
\code{sqrt(.Machine$double.eps)} (which is
\Sexpr{format(sqrt(.Machine$double.eps),digits=3)}).  This is a somewhat
conservative tolerance, in many cases I'd guess
\code{options(lfe.eps=1e-4)} may be sufficient.  This will speed up the
centering.

The package is threaded, that is, it may use more than one cpu.  The
number of threads is fetched upon loading the package, from the
environment variable \env{LFE_THREADS} (or \env{OMP_NUM_THREADS}) and
stored by \code{options(lfe.threads=n)}.  This option may be changed prior to
calling \code{\link{felm}}, if so desired.

Threading is only done for the centering; the extraction of the group
effects is not threaded, but it uses any threading in the underlying
blas-library (which is usually controlled by the \env{OMP_NUM_THREADS}
environment variable).

The package has been tested on datasets with approx 20,000,000
observations with approx 2,300,000 and 270,000 group levels (the
\code{\link{felm}} takes 1-2 hours on 8 cpus, the \code{\link{getfe}} takes a
couple of days).  It uses the sparse Cholesky solver of package
\pkg{Matrix}, which relies heavily on the blas-library.
It's thus strongly recommended to link an optimized blas into R (such as
'goto', 'atlas', 'acml' or 'mkl').

The package will work with any positive number of grouping factors, but if
more than two, their interpretation is in general not well understood.
}

\references{
Abowd, J.M. and Kramarz, F. and Margolis, D.N. (1999) \cite{High Wage Workers and High
Wage Firms}, Econometrica 67 (1999), no. 2, 251--333.
  
Abowd, J. and Creecy, R. and Kramarz, F. (2002) \cite{Computing Person
  and Firm Effects Using Linked Longitudinal Employer-Employee
  Data.} Technical Report TP-2002-06, U.S. Census Bureau.
  \url{http://lehd.did.census.gov/led/library/techpapers/tp-2002-06.pdf}

Andrews, M. and Gill, L. and Schank, T. and Upward, R. (2008)
\cite{High wage workers and low wage firms: negative assortative
  matching or limited mobility bias?}
  J.R. Stat. Soc.(A) 171(3), 673--697. 
  \url{http://dx.doi.org/10.1111/j.1467-985X.2007.00533.x}

Gaure, S. (2011) \cite{OLS with Multiple High Dimensional Category Dummies}, Gaure, S. (to appear)
}

\examples{
  x <- rnorm(100)
  x2 <- rnorm(length(x))
  id <- factor(sample(10,length(x),replace=TRUE))
  firm <- factor(sample(3,length(x),replace=TRUE,prob=c(2,1,1)))
  id.eff <- rnorm(nlevels(id))
  firm.eff <- rnorm(nlevels(firm))
  y <- x + 0.25*x2 + id.eff[id] + firm.eff[firm] + rnorm(length(x))
  dset <- data.frame(y,x,x2,id,firm)
  est <- felm(y ~ x+x2,fl=list(id=id,firm=firm),data=dset)
  print(est)
  print(getfe(est))
# compare with an ordinary lm
  summary(lm(y ~ x+x2+id+firm-1,data=dset))
}
\keyword{regression}
\keyword{models}
