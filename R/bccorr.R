orthonormalize <- function(V) {
  structure(V %*% solve(chol(crossprod(V))), ortho=TRUE)
}

narowsum <- function(x, group) {
  opt <- options(warn=-1)
  res <- try(rowsum(x,group), silent=TRUE)
  options(opt)
  narow <- is.na(rownames(res))
  if(any(narow))
      res[!narow,,drop=FALSE]
  else
      res
}
# try without reference
narowsum <- rowsum
# to estimate the bias of the D-variance we need to find tr(D' P_1 D (D' P_{X,F} D)^{-1})
# where P_1 projects out the mean, P_{X,F} projects out X and F. Let's make function
# which evaluates the matrix inside the trace on a vector.

# Now, note that for matrix x, D %*% x is the same as x[f1,]
# whereas t(D) %*% x is the same as rowsum(x,f1). Similarly with F and f2.

# bias corrected variances

bccorr <- function(est, alpha=getfe(est), corrfactors=1L:2L,
                   nocovar=is.null(est$X) && length(est$fe)==2,
                   tol=0.01, maxsamples=Inf) {

  if(length(est$fe) == 2 && nlevels(est$cfactor) != 1) stop('Bias correction only makes sense on data with 1 component')
  if(length(est$fe) < 2) stop('Bias correction only makes sense for two factors')
  if(length(corrfactors) != 2) stop('corrfactors must have length 2')
  if(is.character(corrfactors)) corrfactors <- match(corrfactors, names(est$fe))
  if(min(corrfactors) < 1 || max(corrfactors) > length(est$fe))
      stop('corrfactors specifies too small or large index')

  fe <- est$fe
  f1 <- fe[[corrfactors[[1]]]]
  f2 <- fe[[corrfactors[[2]]]]
  nf1 <- names(fe)[[corrfactors[[1]]]]
  nf2 <- names(fe)[[corrfactors[[2]]]]
  if(!is.null(attr(f1,'x'))) stop('Interacted factors "',nf1,'" are not supported')
  if(!is.null(attr(f2,'x'))) stop('Interacted factors "',nf2,'" are not supported')
  if(!all(c('fe','effect') %in% colnames(alpha)))
      stop('alpha must contain columns "effect" and "fe"')

  d1 <- alpha[alpha['fe']==nf1,'effect'][f1]
  d2 <- alpha[alpha['fe']==nf2,'effect'][f2]

  var1 <- var(d1)
  var2 <- var(d2)
  cov12 <- cov(d1,d2)

  biascorr <- bcvar(if(nocovar && length(fe) > 2) fe[corrfactors] else fe,
                    if(nocovar) NULL else est$X,
                    corrfactors,
                    var1, var2, cov12,
                    sum(est$residuals^2)/est$df,
                    tol,maxsamples)
  delta1 <- biascorr[['delta1']]
  delta2 <- biascorr[['delta2']]
  delta12 <- biascorr[['delta12']]

  vartheta <- var1 - delta1
  varpsi <- var2 - delta2
  covtp <- cov12 - delta12

  c(corr=covtp/sqrt(vartheta*varpsi),
    v1=vartheta,
    v2=varpsi,
    cov=covtp,
    d1=delta1,
    d2=delta2,
    d12=delta12)

}

varbias <- function(index,fe,X,tol=0.01,epsvar,bvar, maxsamples=Inf) {
  if(length(index) != 1) stop("index must have length 1")
  f <- fe[[index]]
  name <- names(fe)[[index]]
  N <- length(f)
  if(!is.null(X) && nrow(X) != N) stop("X must be of same length as factor")
  nlev <- nlevels(f)
  restf <- fe[-index]

# First, make a factor list for projecting out things.
  fmean <- factor(rep(1,N))
  lmean <- list(fmean)

# project out X and F, i.e. everything but f1
# use M_{F,X} = M_F M_{M_F X}
# precompute and orthonormalize M_F X

  if(!is.null(X)) {
    MFX <- list(structure(fmean,x=orthonormalize(demeanlist(X,restf))))
    invfun <- function(v) {
      rowsum(demeanlist(demeanlist(v[f,],MFX), restf), f)
    }
  } else {
    invfun <- function(v) {
      rowsum(demeanlist(v[f,], restf), f)
    }
  }
  

  # We want precision in the final estimate to be, say, 1%
  # Since we subtract the bias from the biased estimate, the precision
  # in the trace computation depends on the current value of the trace
  # i.e. absolute precision should be 0.01*(sum(d1^2) - epsvar*tr)/epsvar

  epsfun <- function(tr) -tol*abs(bvar*N - epsvar*tr)/epsvar
  # the tolerance before mctrace has got a clue about where we are
  # is a problem. If the bias is very large compared to the variance, we will
  # be in trouble. 
  trfun <- function(x,trtol) {
    # return crude estimate of the trace
    if(trtol == 0) return(abs(nlev))
    DtM1x <- rowsum(demeanlist(x,lmean), f)
    # we use absolute tolerance, mctrace wil give us a trtol.
    # we divide by the L2-norm of DtM1x, since we take the
    # inner product with this afterwards
    tol1 <- -trtol/sqrt(colSums(DtM1x^2))/2
    v <- cgsolve(invfun, DtM1x, eps=tol1)
    colSums(DtM1x * v)
  }
  attr(trfun,'IP') <- TRUE
  res <- epsvar*mctrace(trfun,N=N,tol=epsfun, trname=paste(name,'bias'),
                 maxsamples=maxsamples)/N
}

# compute variance of the biased variance estimate
# using the formula for the variance of a quadratic form with normal distribution
# var(x^t A x) = 2 tr(AVAV) + 4mu^t*AVA*mu
# where V is the variance matrix of x, assumed to be sigma^2 I, and mu is the
# expectation of x (i.e. Dtheta).
varvar <- function(index, fe, X, pointest, resvar, tol=0.01) {
  f <- fe[[index]]
  N <- length(f)
  fmean <- list(factor(rep(1,N)))

  if(is.null(X)) {
      MFX <- fe[-index]
      invfun <- function(x) {
        rowsum(demeanlist(x[f,], MFX),f)
      }
  } else {
#    M_{F,X} = M_F M_{M_F X}
    restf <- fe[-index]
    MFX <- list(structure(factor(rep(1,N)), x=orthonormalize(demeanlist(X,restf))))
    invfun <- function(x) {
      rowsum(demeanlist(demeanlist(x[f,],MFX), restf), f)
    }
  }

  Dtheta <- pointest[f]
  DtM1D <- rowsum(demeanlist(Dtheta,fmean), f)
  v <- cgsolve(invfun, DtM1D, eps=tol/4/resvar)
  meanpart <- 4*resvar * sum(DtM1D * v)
  # this mean part is biased upwards. We should correct it.
  # it turns out that we can do this by changing the sign of the
  # trace term, the bias is the same expression as the trace part
  message('mean part=',meanpart/N^2)
  mytol <- -meanpart/10
  trfun <- function(x,trtol) {
    v <- demeanlist(cgsolve(invfun, rowsum(demeanlist(x,fmean),f),
                            eps=mytol/resvar^2/2)[f,],fmean)
    colSums(v * v)
  }
  attr(trfun,'IP') = TRUE
  trpart <- 2*resvar^2 * mctrace(trfun, N=length(f))
  message('mean part=', meanpart/N^2, ' trpart=',trpart/N^2)
  (-trpart+meanpart)/N^2
}

# if positive tolerance, the tolerance is relative to the bias corrected
# covariance. In this case, the biased covariance (bcov) and residual
# variance (epsvar) must be specified. A negative tolerance is an
# absolute tolerance
covbias <- function(index,fe,X,tol=-0.01, epsvar, bcov, maxsamples=Inf) {
  if(length(index) != 2) stop("index must have length 2")
  if(length(fe) < 2) stop("fe must have length >= 2")
  if(tol > 0 && missing(bcov))
      stop("with relative tolerance, bcov must be specified")

  f1 <- fe[[index[[1]]]]
  f2 <- fe[[index[[2]]]]
  nlev1 <- nlevels(f1)
  nlev2 <- nlevels(f2)
  N <- length(f1)
  name <- paste('cov(',paste(names(fe)[index],collapse=','),')',sep='')
  no2list <- fe[-index[[2]]]
  restf <- fe[-index]
  fmean <- factor(rep(1,N))
  lmean <- list(fmean)

  if(length(f2) != N) stop("error factors must have the same length")
  if(!is.null(X)) {
    if(nrow(X) != N) stop("X must have same length as factors")
    MDX <- list(structure(fmean,x=orthonormalize(demeanlist(X,no2list))))
    invfun <- function(v) {
      rowsum(demeanlist(demeanlist(v[f2,],MDX), no2list), f2)
    }
    if(length(restf) > 0) {
      MX <- list(structure(fmean,x=orthonormalize(demeanlist(X,restf))))
      MXfun <- function(v) demeanlist(demeanlist(v, MX), restf)
      invfunX <- function(v) {
        rowsum(MXfun(v[f1,]), f1)
      }
    } else {
      MX <- list(structure(fmean, x=orthonormalize(X)))
      MXfun <- function(v) demeanlist(v, MX)
      invfunX <- function(v) {
        rowsum(MXfun(v[f1,]), f1)
      }
    }
  } else {
    invfun <- function(v) {
      rowsum(demeanlist(v[f2,], no2list), f2)
    }
    MXfun <- function(v) demeanlist(v, restf)
    invfunX <- function(v) {
      rowsum(MXfun(v[f1,]), f1)
    }
  }
  
  trfun <- function(x,trtol) {
    # return crude estimate of the trace
    if(trtol == 0) return(-abs(nlev1-nlev2))
    M1x <- demeanlist(x,lmean)
    DtM1x <- rowsum(M1x,f1)
    FtM1x <- rowsum(M1x,f2)
    d1 <- sqrt(colSums(DtM1x^2))
    d2 <- sqrt(colSums(FtM1x^2))
#    tol12 <-  -(if(trtol==0) abs(tol)*(nlev1+nlev2) else trtol)
    tol12 <- -trtol
    v <- cgsolve(invfunX, DtM1x, eps=tol12/(d1+d2)/3)
    MXv <- rowsum(MXfun(v[f1,]), f2)
    w <- cgsolve(invfun, FtM1x, eps=tol12/sqrt(colSums(MXv^2))/3)
    -colSums(w* MXv)
  }
  # our function does the inner product, not just matrix application. Signal to mctrace.
  attr(trfun,'IP') <- TRUE
  
  if(tol < 0) {
    # scale by N and epsvar since it's trace level precision
    eps <- tol*N/epsvar
  } else {
    # relative precision
    eps <- function(tr) -tol*abs(bcov*N - epsvar*tr)/epsvar
  }
  epsvar*mctrace(trfun, N=N, tol=eps, trname=name,
                 maxsamples=maxsamples)/N
}

bcvar <- function(fe, X, corrfactors, 
                  var1, var2, cov12, epsvar, tol=0.01, maxsamples=Inf) {
  message('1st variance')
  delta1 <- varbias(corrfactors[[1]],fe,X,tol,epsvar,var1,maxsamples)
  if(is.null(X) && length(fe) == 2) {
    f1 <- fe[[corrfactors[[1]]]]
    f2 <- fe[[corrfactors[[2]]]]
    N <- length(f1)
    delta2 <- delta1 - epsvar*(nlevels(f1)-nlevels(f2))/N
    delta12 <- epsvar*nlevels(f1)/N - delta1
  } else {
    message('2nd variance')
    delta2 <- varbias(corrfactors[[2]],fe,X,tol,epsvar,var2,maxsamples)
    eps <- -tol*sqrt((var1-delta1)*(var2-delta2))
    message('covariance')
    delta12 <- covbias(corrfactors,fe,X,eps,epsvar,maxsamples=maxsamples)
  }
  c(var1=var1-delta1, var2=var2-delta2, cov12=cov12-delta12,
    delta1=delta1, delta2=delta2, delta12=delta12)
}


