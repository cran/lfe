\name{print.lm.fe}
\alias{print.lm.fe}
\title{Printing result}
\description{
\code{print} method for class \code{"fe.lm"}, typically objects returned
by \code{\link{felm}}.
}
\usage{
  \method{print}{lm.fe}(x,digits=max(3,getOption('digits')-3),...)
}

\arguments{
  \item{x}{an object of class '"fe.lm"', typically the result of a call
  to \link{felm}.}
  \item{digits}{minimal number of _significant_ digits.}
  \item{...}{other arguments.}
}
\seealso{\code{\link{getfe}}}
