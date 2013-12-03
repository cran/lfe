\name{demeanlist}
\alias{demeanlist}
\title{Centre vectors on multiple groups}
\description{
  Uses the method of alternating projections to centre
  a (model) matrix on multiple groups, as specified by a list of factors.
  This function is called by \code{\link{felm}}, but it has been
  made available as standalone in case it's needed.
}

\usage{
demeanlist(mtx, fl, icpt=0, eps=getOption('lfe.eps'),
           threads=getOption('lfe.threads'),
           progress=getOption('lfe.pint'),
           accel=getOption('lfe.accel'),
           randfact=TRUE)
}

\arguments{
\item{mtx}{matrix whose columns form vectors to be group-centred. mtx
  may also be a list of vectors or matrices.}
\item{fl}{list of factors defining the grouping structure}
\item{icpt}{the position of the intercept, this column is removed from
  the result matrix}
\item{eps}{a tolerance for the centering}
\item{threads}{an integer specifying the number of threads to use}
\item{progress}{integer. If positive, make progress reports (whenever a
  vector is centered, but not more often than every \code{progress} minutes)}
\item{accel}{integer. Set to 1 if Gearhart-Koshy acceleration should be done.}
\item{randfact}{logical. Should the order of the factors be randomized?
  This may improve convergence.}
}

\details{
For each column \code{y} in \code{mtx}, the equivalent of the
following centering is performed, with \code{cy} as the result.
\preformatted{  
cy <- y; oldy <- y-1
while(sqrt(sum((cy-oldy)**2)) >= eps) {
  oldy <- cy
  for(f in fl) cy <- cy - ave(cy,f)
}
}

Beginning with version 1.6, each factor in \code{fl} may contain an
attribute \code{'x'} which is a numeric vector of the same length as
the factor. The centering is then not done on the means of each group,
but on the projection onto the covariate in each group.  That is, with a
covariate \code{x} and a factor \code{f}, it is like projecting out the
interaction \code{x:f}.
}


\value{
If \code{mtx} is a matrix, a matrix of the same shape, possibly with 
column \code{icpt} deleted.
If \code{mtx} is a list of vectors and matrices, a list of the same
length is returned, with the same vector and matrix-pattern, but the
matrices have the column \code{icpt} deleted.
}
\note{
In the case that the design-matrix is too large for R, i.e. with more
than 2 billion entries, it is possible to create a list of
column-vectors instead (provided the vector-length is smaller than 2
billion).  \code{demeanlist} will be able to centre these vectors.

The \code{accel} argument enables Gearhart-Koshy acceleration as
described in theorem 3.16 by Bauschke, Deutsch, Hundal and Park in "Accelerating the
convergence of the method of alternating projections",
Trans. Amer. Math. Soc. 355 pp 3433-3461 (2003).
}

\examples{
## create a 15x3 matrix
mtx <- matrix(rnorm(45),15,3)

## a list of factors
fl <- list(g1=factor(sample(2,nrow(mtx),replace=TRUE)),
           g2=factor(sample(3,nrow(mtx),replace=TRUE)))

## centre on both means and print result
mtx0 <- demeanlist(mtx,fl)
cbind(mtx0,g1=fl[[1]],g2=fl[[2]],comp=compfactor(fl))

for(i in 1:ncol(mtx0))
   for(n in names(fl))
    cat('col',i,'group',n,'level sums:',tapply(mtx0[,i],fl[[n]],mean),'\n')

}