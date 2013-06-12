\name{kaczmarz}
\alias{kaczmarz}
\title{Solve a linear system defined by factors}
\description{
  Uses the Kaczmarz method to solve a system of the type Dx = R, where
  D is the matrix of dummies created from a list of factors.  
}

\usage{
kaczmarz(fl,R,eps=getOption('lfe.eps'),
         init=NULL,threads=getOption('lfe.threads'))
}

\arguments{
\item{fl}{A list of arbitrary factors of the same length}
\item{R}{numeric.  A vector, matrix or list of such of the same length as the factors}
\item{eps}{a tolerance for the method}
\item{init}{numeric. A vector to use as initial value for the Kaczmarz
  iterations. The algorithm converges to the solution closest to this}
\item{threads}{integer. The number of threads to use when \code{R} is
  more than one vector}
}

\value{
A vector \code{x} of length equal to the sum of the number of levels of
the factors in \code{fl}, which solves the system Dx=R. If the system
is inconsistent, the algorithm may not converge, it will give a warning
and return something which may or may not be close to a solution.
By setting \code{eps=0}, maximum accuracy (with convergence warning)
will be achieved.
}

\note{
This function is used by \code{\link{getfe}}, it's quite specialized,
but it might be useful for other purposes too.  

In case of convergence problems, setting \code{options(lfe.usecg=TRUE)} will
cause the kaczmarz() function to dispatch to the conjugate gradient method
of \pkg{Rcgmin}.  This may or may not be faster.
}

\examples{
## create factors
  f1 <- factor(sample(24000,100000,replace=TRUE))
  f2 <- factor(sample(20000,length(f1),replace=TRUE))
  f3 <- factor(sample(10000,length(f1),replace=TRUE))
  f4 <- factor(sample(8000,length(f1),replace=TRUE))
## the matrix of dummies
  D <- t(rBind(as(f1,'sparseMatrix'),as(f2,'sparseMatrix'),
as(f3,'sparseMatrix'), as(f4,'sparseMatrix')))
  dim(D)
## an x
  truex <- runif(ncol(D))
## and the right hand side
  R <- as.vector(D \%*\% truex)
## solve it
  sol <- kaczmarz(list(f1,f2,f3,f4),R)
## verify that the solution solves the system Dx = R
  res <- D \%*\% sol - R
  sqrt(sum(res^2))
## but the solution is not equal to the true x, because the system is
## underdetermined
  sqrt(sum((sol - truex)^2))
## moreover, the solution from kaczmarz has smaller norm
  sqrt(sum(sol^2))
  sqrt(sum(truex^2))
}

