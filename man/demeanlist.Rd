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
demeanlist(mtx,fl,icpt=0,eps=getOption('lfe.eps'),threads=getOption('lfe.threads'))
}

\arguments{
\item{mtx}{matrix whose columns form vectors to be group-centred. mtx
  may also be a list of vectors.}
\item{fl}{list of factors defining the grouping structure}
\item{icpt}{the position of the intercept, this column is removed from the result}
\item{eps}{a tolerance for the centering}
\item{threads}{an integer specifying the number of threads to use}
}

\details{
For each column \code{y} in \code{mtx}, the equivalent of the
following centering is performed, with \code{cy} as the result.
\preformatted{  
cy <- y; oldy <- y-1
while(sqrt(sum((cy-oldy)**2)) >= eps) {
  oldy <- cy
  for(f in fl) cy <- cy - tapply(cy,f,mean)[f]
}
}
}


\value{
If \code{mtx} is a matrix, a matrix of the same shape, possibly with the
column \code{icpt} deleted.
If \code{mtx} is a list of vectors, a list of the same length is returned.
}

\examples{
## create a 15x3 matrix
mtx <- matrix(rnorm(45),15,3)

## a list of factors
fl <- list(g1=factor(sample(2,nrow(mtx),replace=TRUE)),
           g2=factor(sample(3,nrow(mtx),replace=TRUE)))

## centre on both means and print result
mtx0 <- demeanlist(mtx,fl)
print(cbind(mtx0,g1=fl[[1]],g2=fl[[2]],comp=compfactor(fl)))

\dontrun{for(i in 1:ncol(mtx0))
   for(n in names(fl))
    cat('col',i,'group',n,'level sums:',tapply(mtx0[,i],fl[[n]],mean),'\n')
}
}