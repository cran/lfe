% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oldfelm.R
\name{..oldfelm}
\alias{..oldfelm}
\title{Fit a linear model with multiple group fixed effects (old interface)}
\usage{
..oldfelm(
  formula,
  data,
  exactDOF = FALSE,
  subset,
  na.action,
  contrasts = NULL,
  ...
)
}
\arguments{
\item{formula}{an object of class '"formula"' (or one that can be coerced to
that class): a symbolic description of the model to be fitted. Similarly to
'lm'.  See Details.}

\item{data}{a data frame containing the variables of the model.}

\item{exactDOF}{logical. If more than two factors, the degrees of freedom
used to scale the covariance matrix (and the standard errors) is normally
estimated. Setting \code{exactDOF=TRUE} causes \code{felm} to attempt to
compute it, but this may fail if there are too many levels in the factors.
\code{exactDOF='rM'} will use the exact method in
\code{Matrix::rankMatrix()}, but this is slower. If neither of these methods
works, it is possible to specify \code{exactDOF='mc'}, which utilizes a
Monte-Carlo method to estimate the expectation E(x' P x) = tr(P), the trace
of a certain projection, a method which may be more accurate than the
default guess.

If the degrees of freedom for some reason are known, they can be specified
like \code{exactDOF=342772}.}

\item{subset}{an optional vector specifying a subset of observations to be
used in the fitting process.}

\item{na.action}{a function which indicates what should happen when the data
contain \code{NA}s.  The default is set by the \code{na.action} setting of
\code{options}, and is \code{na.fail} if that is unset.  The 'factory-fresh'
default is \code{na.omit}.  Another possible value is \code{NULL}, no
action. \code{na.exclude} is currently not supported.}

\item{contrasts}{an optional list. See the \code{contrasts.arg} of
\code{model.matrix.default}.}

\item{...}{other arguments.  \itemize{

\item \code{cmethod} character. Which clustering method to use. Known
arguments are \code{'cgm'} (the default), \code{'cgm2'} (or \code{'reghdfe'},
its alias). These alternate methods will generally
yield equivalent results, except in the case of multiway clustering with few
clusters along at least one dimension.

\item \code{keepX} logical. To include a copy of the expanded data matrix in
the return value, as needed by \code{\link[=bccorr]{bccorr()}} and \code{\link[=fevcov]{fevcov()}}
for proper limited mobility bias correction.

\item \code{keepCX} logical. Keep a copy of the centred expanded data matrix
in the return value. As list elements \code{cX} for the explanatory
variables, and \code{cY} for the outcome.

\item \code{keepModel} logical. Keep a copy of the model frame.

\item \code{nostats} logical. Don't include covariance matrices in the
output, just the estimated coefficients and various descriptive information.
For IV, \code{nostats} can be a logical vector of length 2, with the last
value being used for the 1st stages.  \item \code{psdef} logical. In case of
multiway clustering, the method of Cameron, Gelbach and Miller may yield a
non-definite variance matrix. Ordinarily this is forced to be semidefinite
by setting negative eigenvalues to zero. Setting \code{psdef=FALSE} will
switch off this adjustment.  Since the variance estimator is asymptotically
correct, this should only have an effect when the clustering factors have
very few levels.

\item \code{kclass} character. For use with instrumental variables. Use a
k-class estimator rather than 2SLS/IV. Currently, the values \verb{'nagar', 'b2sls', 'mb2sls', 'liml'} are accepted, where the names are from
\cite{Kolesar et al (2014)}, as well as a numeric value for the 'k' in
k-class. With \code{kclass='liml'}, \code{felm} also accepts the argument
\verb{fuller=<numeric>}, for using a Fuller adjustment of the
liml-estimator.

\item \verb{Nboot, bootexpr, bootcluster} Since \code{felm} has quite a bit
of overhead in the creation of the model matrix, if one wants confidence
intervals for some function of the estimated parameters, it is possible to
bootstrap internally in \code{felm}.  That is, the model matrix is resampled
\code{Nboot} times and estimated, and the \code{bootexpr} is evaluated
inside an \code{sapply}. The estimated coefficients and the left hand
side(s) are available by name. Any right hand side variable \code{x} is
available by the name \code{var.x}.  The \code{"felm"}-object for each
estimation is available as \code{est}. If a \code{bootcluster} is specified
as a factor, entire levels are resampled. \code{bootcluster} can also be a
function with no arguments, it should return a vector of integers, the rows
to use in the sample. It can also be the string 'model', in which case the
cluster is taken from the model. \code{bootexpr} should be an expression,
e.g. like \code{quote(x/x2 * abs(x3)/mean(y))}.  It could be wise to specify
\code{nostats=TRUE} when bootstrapping, unless the covariance matrices are
needed in the bootstrap. If you need the covariance matrices in the full
estimate, but not in the bootstrap, you can specify it in an attribute
\code{"boot"} as \code{nostats=structure(FALSE, boot=TRUE)}.

\item \verb{iv, clustervar} deprecated.  These arguments will be removed at
a later time, but are still supported in this field. Users are
\emph{STRONGLY} encouraged to use multipart formulas instead.  In
particular, not all functionality is supported with the deprecated syntax;
iv-estimations actually run a lot faster if multipart formulas are used, due
to new algorithms which I didn't bother to shoehorn in place for the
deprecated syntax.

}}
}
\description{
Fit a linear model with multiple group fixed effects (old interface)
}
