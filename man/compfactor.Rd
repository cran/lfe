\name{compfactor}
\alias{compfactor}
\title{Finding the connection components}
\description{
'compfactor' computes the connection components of the dummy-part of the
model with multiple fixed effects.
}
\usage{
  compfactor(fl)
}

\arguments{
  \item{fl}{a list of factors defining the fixed effects}
}
\value{
A factor of the same length as the factors in the input argument.  It
defines the connection components.
E.g. \code{nlevels(compfactor(fl))} will yield the number of
connection components.
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