#' Overview. Linear Group Fixed Effects
#'
#' The package uses the Method of Alternating Projections to estimate linear
#' models with multiple group fixed effects.  A generalization of the within
#' estimator. It supports IV-estimation with multiple endogenous variables via
#' 2SLS, with conditional F statistics for detection of weak instruments. It is
#' thread-parallelized and intended for large problems. A method for correcting
#' limited mobility bias is also included.
#'
#'
#' This package is intended for linear models with multiple group fixed
#' effects, i.e. with 2 or more factors with a large number of levels.  It
#' performs similar functions as [stats::lm()], but it uses a special
#' method for projecting out multiple group fixed effects from the normal
#' equations, hence it is faster. It is a generalization of the within
#' estimator.  This may be required if the groups have high cardinality (many
#' levels), resulting in tens or hundreds of thousands of dummy variables.  It
#' is also useful if one only wants to control for the group effects, without
#' actually estimating them.  The package may optionally compute standard
#' errors for the group effects by bootstrapping, but this is a very time- and
#' memory-consuming process compared to finding the point estimates.  If you
#' only have a single huge factor, the package \pkg{plm} is probably better
#' suited.  If your factors don't have thousands of levels,
#' [stats::lm()] or other packages are probably better suited.
#' \pkg{lfe} is designed to produce the same results as [stats::lm()]
#' will do if run with the full set of dummies.
#'
#' Projecting out interactions between continuous covariates and factors is
#' supported. I.e. individual slopes, not only individual intercepts. Multiple
#' left hand sides are supported.
#'
#' The package does not support non-linear models. For GLMs with many dummies there
#' is a package \pkg{alpaca} which uses similar methods to project them out.
#'
#' The estimation is done in two steps.  First the other coefficients are
#' estimated with the function [felm()] by centering on all the group
#' means, followed by an OLS (similar to lm).  Then the group effects are
#' extracted (if needed) with the function [getfe()].  This method is
#' described by \cite{Gaure (2013)}, but also appears in \cite{Guimaraes and
#' Portugal (2010)}, disguised as the Gauss-Seidel algorithm.
#'
#' There's also a function [demeanlist()] which just does the
#' centering on an arbitrary matrix or data frame, and there's a function
#' [compfactor()] which computes the connected components which are
#' used for interpreting the group effects when there are only two factors (see
#' the Abowd et al references), they are also returned by [getfe()].
#'
#' For those who study the correlation between the fixed effects, like in
#' \cite{Abowd et al. (1999)}, there are functions [bccorr()] and
#' [fevcov()] for computing limited mobility bias corrected
#' correlations and variances with the method described in \cite{Gaure
#' (2014b)}.
#'
#' Instrumental variable estimations are supported with 2SLS. Conditional F
#' statistics for testing reduced rank weak instruments as in \cite{Sanderson
#' and Windmeijer (2015)} are available in [condfstat()].  Joint
#' significance testing of coefficients is available in [waldtest()].
#'
#' The centering on the means is done with a tolerance which is set by
#' `options(lfe.eps=1e-8)` (the default).  This is a somewhat conservative
#' tolerance, in many cases I'd guess `1e-6` may be sufficient.  This may
#' speed up the centering.  In the other direction, setting
#' `options(lfe.eps=0)` will provide maximum accuracy at the cost of
#' computing time and warnings about convergence failure.
#'
#' The package is threaded, that is, it may use more than one cpu.  The number
#' of threads is fetched upon loading the package from the environment variable
#' \env{LFE_THREADS}, \env{OMP_THREAD_LIMIT}, \env{OMP_NUM_THREADS} or
#' \env{NUMBER_OF_PROCESSORS} (for Windows), and stored by
#' `options(lfe.threads=n)`.  This option can be changed prior to calling
#' [felm()], if so desired.  Note that, typically, \pkg{lfe} is
#' limited by memory bandwidth, not cpu speed, thus fast memory and large cache
#' is more important than clock frequency. It is therefore also not always true
#' that running on all available cores is much better than running on half of
#' them.
#'
#' Threading is only done for the centering; the extraction of the group
#' effects is not threaded. The default method for extracting the group
#' coefficients is the iterative Kaczmarz-method, its tolerance is also the
#' `lfe.eps` option. For some datasets the Kaczmarz-method is converging
#' very slowly, in this case it may be replaced with a conjugate gradient
#' method by setting the option `options(lfe.usecg=TRUE)`. Various
#' time-consuming parts of \pkg{lfe} may print progress reports, the minimum
#' interval in seconds is `options(lfe.pint=1800)`.
#'
#' The package has been tested on datasets with approx 20,000,000 observations
#' with 15 covariates and approx 2,300,000 and 270,000 group levels (the
#' [felm()] took about 50 minutes on 8 cpus, the [getfe()]
#' takes 5 minutes).  Though, beware that not only the size of the dataset
#' matters, but also its structure, as demonstrated by \cite{Gaure (2014a)}.
#'
#' The package will work with any number of grouping factors, but if more than
#' two, their interpretation is in general not well understood, i.e. one should
#' make sure that the group coefficients are estimable. A discussion of
#' estimability, the algorithm used, and convergence rate are available in
#' vignettes, as well as in the published papers in the citation list
#' (`citation('lfe')`).
#'
#' In the exec-directory there is a perl-script `lfescript` which is used
#' at the author's site for automated creation of R-scripts from a simple
#' specification file.  The format is documented in `doc/lfeguide.txt`.
#'
#' \pkg{lfe} is similar in function, though not in method, to the Stata modules
#' `a2reg` and `felsdvreg`.  The method is very similar to the one in
#' the Stata module `reghdfe`.
#'
#' @name lfe-package
#' @aliases lfe-package lfe
#' @references Abowd, J.M., F. Kramarz and D.N. Margolis (1999) \cite{High Wage
#' Workers and High Wage Firms}, Econometrica 67 (1999), no. 2, 251--333.
#' \doi{10.1111/1468-0262.00020}
#'
#' Abowd, J.M., R. Creecy and F. Kramarz (2002) \cite{Computing Person and Firm
#' Effects Using Linked Longitudinal Employer-Employee Data.} Technical Report
#' TP-2002-06, U.S. Census Bureau.
#' <https://www2.census.gov/ces/tp/tp-2002-06.pdf>
#'
#' Andrews, M., L. Gill, T. Schank and R. Upward (2008) \cite{High wage workers
#' and low wage firms: negative assortative matching or limited mobility bias?}
#' J.R. Stat. Soc.(A) 171(3), 673--697.
#' \doi{10.1111/j.1467-985X.2007.00533.x}
#'
#' Cornelissen, T. (2008) \cite{The stata command felsdvreg to fit a linear
#' model with two high-dimensional fixed effects.} Stata Journal,
#' 8(2):170--189, 2008.
#' <https://econpapers.repec.org/RePEc:tsj:stataj:v:8:y:2008:i:2:p:170-189>
#'
#' Correia, S. (2014) \cite{REGHDFE: Stata module to perform linear or
#' instrumental-variable regression absorbing any number of high-dimensional
#' fixed effects}, Statistical Software Components, Boston College Department
#' of Economics. <https://econpapers.repec.org/RePEc:boc:bocode:s457874>
#'
#' Croissant, Y. and G. Millo (2008) \cite{Panel Data Econometrics in R: The
#' plm Package}, Journal of Statistical Software, 27(2).
#' <https://www.jstatsoft.org/v27/i02/>
#'
#' Gaure, S. (2013) \cite{OLS with Multiple High Dimensional Category
#' Variables.} Computational Statistics and Data Analysis, 66:8--18, 2013
#' \doi{10.1016/j.csda.2013.03.024}
#'
#' Gaure, S. (2014a) \cite{lfe: Linear Group Fixed Effects.} The R Journal,
#' 5(2):104-117, Dec 2013.
#' <https://journal.r-project.org/archive/2013/RJ-2013-031/RJ-2013-031.pdf>
#'
#' Gaure, S. (2014b), \cite{Correlation bias correction in two-way
#' fixed-effects linear regression}, Stat 3(1):379-390, 2014.
#' \doi{10.1002/sta4.68}
#'
#' Guimaraes, P. and Portugal, P. (2010) \cite{A simple feasible procedure to
#' fit models with high-dimensional fixed effects.} The Stata Journal,
#' 10(4):629--649, 2010.
#' <https://www.stata-journal.com/article.html?article=st0212>
#'
#' Ouazad, A. (2008) \cite{A2REG: Stata module to estimate models with two
#' fixed effects.} Statistical Software Components S456942, Boston College
#' Department of Economics.
#' <https://ideas.repec.org/c/boc/bocode/s456942.html>
#'
#' Sanderson, E. and F. Windmeijer (2014) \cite{A weak instrument F-test in
#' linear IV models with multiple endogenous variables}, Journal of
#' Econometrics, 2015.
#' <https://www.sciencedirect.com/science/article/pii/S0304407615001736>
#' @keywords regression models
#' @examples
#'
#' oldopts <- options("lfe.threads")
#' options(lfe.threads = 2)
#' x <- rnorm(1000)
#' x2 <- rnorm(length(x))
#' id <- factor(sample(10, length(x), replace = TRUE))
#' firm <- factor(sample(3, length(x), replace = TRUE, prob = c(2, 1.5, 1)))
#' year <- factor(sample(10, length(x), replace = TRUE, prob = c(2, 1.5, rep(1, 8))))
#' id.eff <- rnorm(nlevels(id))
#' firm.eff <- rnorm(nlevels(firm))
#' year.eff <- rnorm(nlevels(year))
#' y <- x + 0.25 * x2 + id.eff[id] + firm.eff[firm] +
#'   year.eff[year] + rnorm(length(x))
#' est <- felm(y ~ x + x2 | id + firm + year)
#' summary(est)
#'
#' getfe(est, se = TRUE)
#' # compare with an ordinary lm
#' summary(lm(y ~ x + x2 + id + firm + year - 1))
#' options(oldopts)
#'
#' @useDynLib lfe, .registration=TRUE, .fixes='C_'
#' @importFrom methods as
#' @importFrom xtable xtable
#' @importFrom sandwich estfun bread
#' @import 'stats'
#' @import Formula
#' @importFrom Matrix t Diagonal rankMatrix Cholesky nnzero crossprod tcrossprod diag
#' @importClassesFrom Matrix sparseMatrix
"_PACKAGE"
