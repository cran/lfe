\name{kaczmarz}
\alias{kaczmarz}
\title{Solve a linear system defined by factors}
\description{
  Uses the Kaczmarz method to solve a system of the type Dx = R, where
  D is the matrix of dummies created from a list of factors.  
}

\usage{
kaczmarz(fl,R,eps=getOption('lfe.eps'))
}

\arguments{
\item{fl}{A list of arbitrary factors of the same length}
\item{R}{numeric.  A vector of the same length as the factors}
\item{eps}{a tolerance for the method}
}

\value{
A vector \code{x} of length equal to the sum of the number of levels of
the factors in \code{fl}, which solves the system Dx=R. If the system
is inconsistent, the algorithm may not converge, it will give a warning
and return something which may or may not be close to a solution.
By setting \code{eps=0}, maximum accuracy (with convergence warning)
will be achieved.

The returned value has a numeric attribute \code{"tolerance"} containing
a measure of the accuracy of the solution.
}

\note{
This function is used by \code{\link{getfe}}, it's quite specialized,
but it might be useful for other purposes too.  
}

\examples{
## create factors
  f1 <- factor(sample(4000,10000,replace=TRUE))
  f2 <- factor(sample(2500,length(f1),replace=TRUE))
  f3 <- factor(sample(1000,length(f1),replace=TRUE))
  f4 <- factor(sample(800,length(f1),replace=TRUE))
## the matrix of dummies
  D <- t(rBind(as(f1,'sparseMatrix'),as(f2,'sparseMatrix'),
as(f3,'sparseMatrix'), as(f4,'sparseMatrix')))
  dim(D)
## an x
  truex <- rnorm(ncol(D))
## and the right hand side
  R <- as.vector(D \%*\% truex)
## solve it
  sol <- kaczmarz(list(f1,f2,f3,f4),R)
## verify that the solution solves the system
  res <- D \%*\% sol - R
  sqrt(sum(res^2))
## but the solution is not equal to the true x, because the system is
## underdetermined
  sqrt(sum((sol - truex)^2))
}

