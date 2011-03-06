\name{compfactor}
\alias{compfactor}
\title{Find the connected components}
\description{
'compfactor' computes the connected components of the dummy-part of the
model.
}
\usage{
  compfactor(fl)
}

\arguments{
  \item{fl}{a list of factors defining the dummies}
}
\value{
A factor of the same length as the factors in the input argument.  It
defines the connected components.
E.g. \code{nlevels(compfactor(fl))} will yield the number of
connected components.
}

\details{
If there are more than two factors, only the first two will be used.
}
\examples{
## create two factors
f1 <- factor(sample(300,400,replace=TRUE))
f2 <- factor(sample(300,400,replace=TRUE))

## find the components
cf <- compfactor(list(f1=f1,f2=f2))

## show the third largest component
fr <- data.frame(f1,f2,cf)
fr[cf==3,]
}