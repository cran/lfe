\name{summary.felm}
\alias{summary.felm}
\title{Summarizing felm model fits}
\description{
\code{summary} method for class \code{"felm"}.
}

\usage{ \method{summary}{felm}(object, ..., robust=!is.null(object$cse)) }

\arguments{
  \item{object}{an object of class \code{"felm"}, a result of a call to \code{felm}.}
  \item{robust}{logical.}
  \item{...}{further arguments passed to or from other methods.}
}

\value{
  The function \code{summary.felm} returns an object of \code{class}
  \code{"summary.felm"}.  
  It is
  quite similar to en \code{"summary.lm"} object, but not entirely compatible.

  The \code{"summary.felm"} object is a list containing the following fields:

  \item{residuals}{a numerical vector. The residuals of the full
  system, with dummies.}
  \item{p}{an integer. The total number of coefficients, including those
    projected out.}
  \item{coefficients}{a p x 4 matrix with columns for the estimated
    coefficients, their standard errors, t-statistic and corresponding
    (two-sided) p-value.}
  \item{rse}{residual standard error.}
  \item{r2}{R^2, the fraction of variance explained by the model.}
  \item{r2adj}{Adjusted R^2.}
  \item{fstat}{F-statistic.}
  \item{pval}{P-values.}
  \item{fe}{list of factors. A list of the \code{G()} terms in the model.}

}

\note{
The standard errors are adjusted for the reduced degrees of freedom
coming from the dummies which are implicitly present.  They are also
small-sample corrected.

If the \code{robust} parameter is \code{FALSE}, the returned object will
contain ordinary standard errors. If the \code{robust} parameter is
\code{TRUE}, clustered standard errors are reported if a
\code{clustervar} was specified in the call to \code{felm}; if not,
heteroskedastic robust standard error are reported.
}

