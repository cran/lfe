\name{summary.felm}
\alias{summary.felm}
\title{Summarize felm model fits}
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
  \item{fe}{list of factors. A list of the terms in the second part of
  the model.}
  \item{iv1fstat}{F-statistic for excluded instruments in 1. step IV, see \code{\link{felm}}.}
  \item{iv1pval}{P-value for \code{iv1fstat}.}

}

\note{
The standard errors are adjusted for the reduced degrees of freedom
coming from the dummies which are implicitly present.  They are also
small-sample corrected.

If the \code{robust} parameter is \code{FALSE}, the returned object will
contain ordinary standard errors. If the \code{robust} parameter is
\code{TRUE}, clustered standard errors are reported if a cluster was specified in the call to \code{felm}; if not,
heteroskedastic robust standard error are reported.
The F-statistic is not reliable for robust and clustered errors, the
reason being that the author has not yet figured out how to compute it
without the full covariance matrix for the dummies that have been
projected out.  In this case, the \code{print} method will report the F-statistic
under the assumption of normally distributed residuals (which can be
computed from the residuals only), and an
F-statistic (via a Wald-statistic) for robust/clustered errors which does not take into account
the covariates which have been projected out.

For a 1st stage IV-regression, an F-statistic against the model with
omitted instruments is also computed, but this one also assumes that the
coefficients for the absorbed dummies are equal, so
that the covariance comes from the ordinary covariates only.
}

