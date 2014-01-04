\name{felm}
\alias{felm}
\title{Fitting linear models with multiple group fixed effects}
\description{
'felm' is used to fit linear models with multiple group fixed effects,
similarly to lm.  It uses the Method of Alternating projections to
sweep out multiple group effects from the normal equations before
estimating the remaining coefficients with OLS.

This function is intended for use with large datasets with multiple group
effects of large cardinality.  If dummy-encoding the group effects
results in a manageable number of coefficients, you are probably better
off by using \code{\link{lm}}.
}

\usage{ felm(formula, data, iv=NULL, clustervar=NULL, exactDOF=FALSE,
subset, na.action, contrasts=NULL) }

\arguments{
  \item{formula}{an object of class '"formula"' (or one that can be
  coerced to that class: a symbolic description of the model to be
  fitted. Similarly to 'lm'.  Grouping factors \code{f} are coded as G(f).
  Interactions between a covariate \code{x} and a factor \code{f} can be
  projected out with the syntax \code{G(x:f)}.
}

  \item{data}{a data frame containing the variables of the model}

  \item{iv}{a formula describing an instrumented variable. Estimated via two step OLS}
  
  \item{clustervar}{a string or factor.  Either the name of a variable
        or a factor. Used for computing clustered standard errors.}

  \item{exactDOF}{logical. If more than two factors, the degrees
        of freedom used to scale the covariance matrix (and the standard
        errors) is normally estimated. Setting \code{exactDOF=TRUE}
        causes \code{felm} to attempt to compute it, but this may fail if there are
        too many levels in the factors.  \code{exactDOF='rM'} will use
        the exact method in Matrix::rankMatrix(), but this is slower. If
        the degrees of freedom for some reason are known, they can be
        specified like \code{exactDOF=342772}} 

 \item{subset}{an optional vector specifying a subset of observations to
   be used in the fitting process.}
 
 \item{na.action}{a function which indicates what should happen when the
        data contain \code{NA}s.  The default is set by the
        \code{na.action} setting of \code{options}, and is
        \code{na.fail} if that is unset.  The 'factory-fresh' default is
        \code{na.omit}.  Another possible value is \code{NULL}, no
        action. \code{na.exclude} is currently not supported.}


  \item{contrasts}{an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}}
}

\value{
  \code{felm} returns an object of \code{class} \code{"felm"}.  It is
  quite similar to an \code{"lm"} object, but not entirely compatible.

  The generic \code{summary}-method will yield a summary which may be \code{print}'ed.
  The object has some resemblance to the an \code{lm} object, and some
  postprocessing methods designed for \code{lm} may happen to work. It
  may however be necessary to coerce the object to succeed with this.

  The \code{"felm"} object is a list containing the following fields:

  \item{coefficients}{a numerical vector. The estimated coefficients.}
  \item{N}{an integer. The number of observations}
  \item{p}{an integer. The total number of coefficients, including those
    projected out.}
  \item{response}{a numerical vector. The response vector.}
  \item{fitted.values}{a numerical vector. The fitted values.}
  \item{residuals}{a numerical vector. The residuals of the full
  system, with dummies.}
  \item{r.residuals}{a numerical vector. Reduced residuals, i.e. the residuals resulting from
  predicting \emph{without} the dummies.}
  \item{cfactor}{factor of length N. The factor describing the connected
    components of the two first \code{G()} terms in the model.}
  \item{vcv}{a matrix. The variance-covariance matrix.}
  \item{fe}{list of factors. A list of the \code{G()} terms in the model.}

  \item{step1}{list of '\code{felm}' objects for the IV 1. step(s), if
    used.}
  
  \item{iv1fstat}{numerical vector. For IV 1. steps, F-value for
    excluded instruments, the number of parameters in restricted model
    and in the unrestricted model.}
    
}

\details{
The \code{G()} is not a function in itself, it is just syntax to
distinguish the grouping factors. It does, however, translate to
\code{as.factor()} inside \code{felm()}.  For the \code{G(x:f)} syntax,
the \code{x} must be numeric vector or matrix or factor, and \code{as.factor()} is
applied to \code{f}.  The entity inside \code{G()} is not treated as an
ordinary formula, in particular it is not possible with things like \code{G(x*f)}.

The standard errors are adjusted for the reduced degrees of freedom
coming from the dummies which are implicitly present.  In the case of
two factors, the exact number of implicit dummies is easy to compute.  If there
are more factors, the number of dummies is estimated by assuming there's
one reference-level for each factor, this may be a slight over-estimation,
leading to slightly too large standard errors. Setting \code{exactDOF}
computes the exact degrees of freedom with \code{rankMatrix()} in package \pkg{Matrix}.
Note that version 1.1-0 og \pkg{Matrix} has a bug in \code{rankMatrix()}
for sparse matrices which may cause cause it to return the wrong value. A fix is underway.

For the \code{iv}-argument, it is only necessary to include the instruments on the
right hand side.  The other covariates, from \code{formula}, are added automatically
in the first step.  See the examples. \code{iv} can also be a list of
formulas if more than one variable is instrumented.  However, all
instruments should then be specified in all the formulas.  A more
consise syntax for multiple instruments will probably be implemented in
the future.

The \code{contrasts} argument is similar to the one in \code{lm()}, it
is used for the factors outside the \code{G()} terms. The factors inside
the \code{G()} terms are analyzed as part of a possible subsequent \code{getfe()} call.

Ideally, the \code{clustervar} should have been an option to the
\code{summary}-function instead.  However, this would require keeping a
copy of the model matrix in the returned structure.  Since this function
is intended for very large datasets, we discard the model matrix to save
memory, keeping only residuals and other summary statistics.

Note that the syntax of the \code{felm}-function has changed, it does no longer
allow a separate specification of the group factors, they must be specified
with the \code{G()}-syntax.  The old felm is still available as \code{lfe:::felm.old}, but
it will no longer be maintained.
}

\seealso{\code{\link{getfe}}}
\examples{
## create covariates
x <- rnorm(1000)
x2 <- rnorm(length(x))

## individual and firm
id <- factor(sample(20,length(x),replace=TRUE))
firm <- factor(sample(13,length(x),replace=TRUE))

## effects for them
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

## left hand side
u <- rnorm(length(x))
y <- x + 0.5*x2 + id.eff[id] + firm.eff[firm] + u

## estimate and print result
est <- felm(y ~ x+x2+G(id)+G(firm))
summary(est)
## compare with lm
summary(lm(y ~ x + x2 + id + firm-1))
## alternatively
\dontrun{felm(y ~ x + x2,fl=list(id=id,firm=firm))
  getfe(est)
}


# make a weird example with 'reverse causation'
# Q and W are instrumented by x3 and G(x4), report robust s.e.
x3 <- rnorm(length(x))
x4 <- sample(12,length(x),replace=TRUE)

Q <- 0.3*x3 + x + 0.2*x2 + id.eff[id] + 0.3*log(x4) - 0.3*y + rnorm(length(x),sd=0.3)
W <- 0.7*x3 - 2*x + 0.1*x2 - 0.7*id.eff[id] + 0.8*cos(x4) - 0.2*y+ rnorm(length(x),sd=0.6)

# add them to the outcome
y <- y + Q + W

ivest <- felm(y ~ x + x2 + G(id)+G(firm) + Q + W, iv=list(Q ~ x3+G(x4), W ~x3+G(x4)))
summary(ivest,robust=TRUE)
\dontrun{
# compare with the not instrumented fit:
summary(felm(y ~ x + x2 + G(id)+G(firm) + Q + W))
}
}
\keyword{regression}
\keyword{models}
