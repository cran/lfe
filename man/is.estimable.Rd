\name{is.estimable}
\alias{is.estimable}
\title{Verify estimability of function}
\description{
  Verify that a function you have written is indeed estimable.
}
\usage{
is.estimable(ef,fe,R=NULL)
}


\arguments{
\item{ef}{function.  The function to be verified.}
\item{fe}{list of factors.}
\item{R}{numeric.  Vector of residuals, is \code{NULL}, one random is created.}
}

\value{
Returns a logical.  
}

\examples{

## create individual and firm
id <- factor(sample(5000,50000,replace=TRUE))
firm <- factor(sample(3000,50000,replace=TRUE))

## create some estimable functions. It's faster to
## use numerical indices in ef rather than strings, and the input v
## to ef has no names, we have to add them when requested
ef <- function(v,addnames) {
  w <- c(v[6]-v[5],v[7000]+v[5],v[7000]-v[6000])
  if(addnames) names(w) <-c('id6-id5','f2k+id5','f2k-f1k')
  w
}
is.estimable(ef,list(id=id,firm=firm))

## Then make an error; in the last coordinate, sum two firms
ef <- function(v,addnames) {
  w <- c(v[6]-v[5],v[7000]+v[5],v[7000]+v[6000])
  if(addnames) names(w) <-c('id6-id5','f2k+id5','f2k-f1k')
  w
}
is.estimable(ef,list(id=id,firm=firm))

}
