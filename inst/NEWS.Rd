\name{NEWS}
\title{lfe news}

\section{Changes in version 2.2 (2015-04-01)}{
  \itemize{

    \item A bug which caused \code{\link{felm}} to crash R when all
    cases of the data had NAs has been fixed. An error is raised instead.

    \item \code{\link{felm}} now supports weights.

    \item The conjugate gradient solver and trace estimation routine
    used by \code{\link{fevcov}} have been made public as
    \code{\link{cgsolve}} and \code{\link{mctrace}}.

    
  }
}

\section{Changes in version 2.0-2.1 (2015-03-05)}{
\itemize{
    \item Better support for matrices as covariates, both on the left
and right hand side of the equation, and in the instrumental variable
specification.

   \item Support for multiple left hand sides, with formulas like
\code{y|z|w ~ X |f1 + f2}.  The content of 'felm' objects have changed
slightly to accomodate this, and generic functions like vcov() and
coef() have gained an extra argument 'lhs' to specify which left hand
side to retrieve information about.

   \item Two new functions, waldtest() and condfstat() have been included to
test linear restrictions on parameters, and to compute conditional F
statistics to test for weak instruments in the IV estimations.
Observations with missing values are now correctly removed in IV-estimations.

   \item The old syntax, with 'G()' to specify factors to project out, has been
deprecated. It still works, but will yield a warning.  Similarly with
the arguments 'iv' and 'clustervar' to felm().

   \item \code{felm} now takes a new argument in the \code{...} part,
\code{kclass} to estimate k-class estimators with instrumental variables.
}
}

\section{Changes in version 1.8 (2014-12-25)}{
\itemize{
  \item Two new functions were added, fevcov() and bccorr(), to compute
limited mobility bias corrections for covariances (and correlations)
between high dimensional factors.  See citation('lfe') for reference
to article describing it.
}
}

\section{Changes in version 1.7 (2014-03-13)}{
\itemize{
   \item  Support for multiway clustering added.  Multiway clusters are specified
in the 4th part of the formula, like \code{y ~ x | f1+f2 | 0 | c1+c2}
}
}

\section{Changes in version 1.6 (2013-12-06)}{
\itemize{
   \item Added support for projecting out interactions between continuous
covariates and factors with specifications like \code{y ~ x | f1 + x2:f1}.

}
}

\section{Version 1.5 (2013-10-28)}{
\itemize{
   \item Version presented in R-journal article. (See citation('lfe') for
reference)
}
}
