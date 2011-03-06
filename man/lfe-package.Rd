\name{lfe-package}
\alias{lfe-package}
\alias{lfe}
\docType{package}
\title{Linear Group Fixed Effects}
\concept{Method of Alternating Projections}
\concept{Kaczmarz Method}
\concept{Fixed Effect Estimator}
\concept{Wihtin Estimator}
\concept{Multiple Fixed Effects}
\description{
The package uses the Method of Alternating Projections to estimate
linear models with multiple group fixed effects.  A generalization
of the within esitmator. It is
thread-parallelized and intended for large problems.
}

\details{

This package is intended for linear models with multiple group fixed effects,
i.e. with 2 or more factors with a large number of levels.
It performs no other functions than \code{\link[stats]{lm}},
but it uses a special method for projecting out multiple group fixed
effects from the normal equations, hence it is faster. It is a
generalization of the within estimator.  This may be required if
the groups have high cardinality (many levels), resulting in tens
or hundreds of thousands of dummy-variables.  It is also useful if one
only wants to control for the group effects, without actually computing
them.  The package may optionally compute standard errors for the group
effects, but this is a very time- and memory-consuming process compared
to finding the point estimates.

The estimation is done in two steps.  First the other coefficients are
estimated with the function \code{\link{felm}} by centering on all the
group means, followed by an OLS (similar to lm).  Then the group effects
are extracted (if needed) with the function \code{\link{getfe}}.  This
method is described in Gaure (2011).

There's also a function \code{\link{demeanlist}} which just does the
centering on an arbitrary matrix, and there's a function
\code{\link{compfactor}} which computes the connected components (which
are used for interpreting the group effects when there are only two
factors, see the Abowd et al references), they are also returned by
\code{\link{getfe}}).

The centering on the means is done with a tolerance which is
set by \code{options(lfe.eps=1e-8)} (the default).  This is a somewhat
conservative tolerance, in many cases I'd guess
\code{1e-6} may be sufficient.  This may speed up the
centering.  In the other direction, setting \code{options(lfe.eps=0)}
will provide maximum accuracy at the cost of computing time and
warnings about convergence failure.

The package is threaded, that is, it may use more than one cpu.  The
number of threads is fetched upon loading the package, from the
environment variable \env{LFE_THREADS} (or \env{OMP_NUM_THREADS}) and
stored by \code{options(lfe.threads=n)}.  This option may be changed prior to
calling \code{\link{felm}}, if so desired.  Note that, typically,
\pkg{lfe} is limited by memory-bandwidth, not cpu-speed, thus fast
memory and large cache is more important than clock-frequency. It's
therefore also not always true that running on all available cores is
much better than running on half of them.

Threading is only done for the centering; the extraction of the group
effects is not threaded. The default method for extracting the
group coefficients is the iterative Kaczmarz-method, its tolerance
is also the \code{lfe.eps} option.

The package has been tested on datasets with approx 20,000,000
observations with 15 covariates and approx 2,300,000 and 270,000 group
levels (the \code{\link{felm}} took about 50 minutes on 8 cpus, the
\code{\link{getfe}} takes 5 minutes).  Though, beware that not
only the size of the dataset matters, but also its structure.

The package will work with any positive number of grouping factors, but if
more than two, their interpretation is in general not well understood,
i.e. one should make sure that the coefficients are estimable.

With certain types of data, which the author has not got a clear
understanding of, the \code{getfe} function with the Kaczmarz method
has convergence problems.  In this case, the option
\code{options(lfe.randkac=TRUE)} may be set prior to calling
\code{getfe}.  This causes random shuffling of the residual system,
which sometimes improves the convergence.

\pkg{lfe} is similar in function, though not in method, to the
Stata modules \code{a2reg} and \code{felsdvreg}.
}

\references{
Abowd, J.M., F. Kramarz and D.N. Margolis (1999) \cite{High Wage Workers and High
Wage Firms}, Econometrica 67 (1999), no. 2, 251--333.
\url{http://dx.doi.org/10.1111/1468-0262.00020}
  
Abowd, J.M., R. Creecy and F. Kramarz (2002) \cite{Computing Person
  and Firm Effects Using Linked Longitudinal Employer-Employee
  Data.} Technical Report TP-2002-06, U.S. Census Bureau.
  \url{http://lehd.did.census.gov/led/library/techpapers/tp-2002-06.pdf}

Andrews, M., L. Gill, T. Schank and R. Upward (2008)
\cite{High wage workers and low wage firms: negative assortative
  matching or limited mobility bias?}
  J.R. Stat. Soc.(A) 171(3), 673--697. 
  \url{http://dx.doi.org/10.1111/j.1467-985X.2007.00533.x}

Cornelissen, T. (2008)
\cite{The stata command felsdvreg to fit a linear model with two
  high-dimensional fixed effects.}
Stata Journal, 8(2):170--189, 2008.
\url{http://econpapers.repec.org/RePEc:tsj:stataj:v:8:y:2008:i:2:p:170-189}

Gaure, S. (2011) \cite{OLS with Multiple High Dimensional Category Variables} (to appear)

Ouazad, A. (2008)
\cite{A2REG: Stata module to estimate models with two fixed effects.}
Statistical Software Components S456942, Boston College Department of Economics.
\url{http://ideas.repec.org/c/boc/bocode/s456942.html}

}

\examples{
  x <- rnorm(1000)
  x2 <- rnorm(length(x))
  id <- factor(sample(10,length(x),replace=TRUE))
  firm <- factor(sample(3,length(x),replace=TRUE,prob=c(2,1.5,1)))
  year <- factor(sample(10,length(x),replace=TRUE,prob=c(2,1.5,rep(1,8))))
  id.eff <- rnorm(nlevels(id))
  firm.eff <- rnorm(nlevels(firm))
  year.eff <- rnorm(nlevels(year))
  y <- x + 0.25*x2 + id.eff[id] + firm.eff[firm] +
         year.eff[year] + rnorm(length(x))
  est <- felm(y ~ x+x2+G(id)+G(firm)+G(year))
  summary(est)

  getfe(est,se=TRUE)
# compare with an ordinary lm
  summary(lm(y ~ x+x2+id+firm+year-1))
}
\keyword{regression}
\keyword{models}
