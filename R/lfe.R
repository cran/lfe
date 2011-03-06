# Author: Simen Gaure
# Copyright: 2011, Simen Gaure
# Licence: Artistic 2.0

.onLoad <- function(libname,pkgname) {
  if(is.null(getOption('lfe.eps')))
    options(lfe.eps=1e-8)
  if(is.null(getOption('lfe.pint')))
    options(lfe.pint=300)
  if(is.null(getOption('lfe.threads'))) {
    cr <- as.integer(Sys.getenv('LFE_THREADS'))
    if(is.na(cr)) cr <- as.integer(Sys.getenv('OMP_NUM_THREADS'))
    if(is.na(cr)) cr <- numcores()
    if(is.na(cr)) {
      packageStartupMessage("LFE running with 1 thread, set options(lfe.threads=n) to use more")
      cr <- 1
    }
    options(lfe.threads=cr)
  }
}


numcores <- function() {
  "This function is snatched from package 'multicore'"
  systems <- list(darwin = "/usr/sbin/sysctl -n hw.ncpu 2>/dev/null",
                  freebsd = "/sbin/sysctl -n hw.ncpu 2>/dev/null",
                  linux = "grep processor /proc/cpuinfo 2>/dev/null|wc -l",
                  irix = c("hinv |grep Processors|sed 's: .*::'", "hinv|grep '^Processor '|wc -l"),
                  solaris = "/usr/sbin/psrinfo -v|grep 'Status of.*processor'|wc -l")
  for (i in seq(systems)) if (length(grep(paste("^",
                                                names(systems)[i], sep = ""), R.version$os)))
    for (cmd in systems[i]) {
      a <- gsub("^ +", "", system(cmd, TRUE)[1])
      if (length(grep("^[1-9]", a)))
        return(as.integer(a))
    }
  NA
}



#   Some things in this file is done in a weird way.
#   In some cases there are efficiency reasons for this, e.g. because
#   the "standard" way of doing things may result in a copy which is costly
#   when the problem is *large*.
#   In other cases it may simply be due to the author's unfamiliarity with how
#   things should be done in R


demeanlist <- function(mtx,fl,icpt=0,eps=getOption('lfe.eps'),
                       threads=getOption('lfe.threads'),
		       progress=getOption('lfe.pint')) {
  if(is.null(threads)) threads <- 1
  islist <- is.list(mtx)
  if(!islist) mtx <- list(mtx)
  res <- .Call('demeanlist',
     mtx,
     as.list(fl),
     as.integer(icpt),               
     as.double(eps),
     as.integer(threads),
     as.integer(progress),
		PACKAGE='lfe')

  if(!islist) res <- res[[1]]
  names(res) <- names(mtx)
  return(res)
}

compfactor <- function(fl) {
  if(length(fl) == 1) return(factor(rep(1,length(fl[[1]]))))
  factor(.Call('conncomp',fl[1:2]))
}


# a function for group fixed effect lm, i.e. first complete centering
# on category variables, then ols

felm <- function(formula,fl,data) {
  mf <- match.call(expand.dots = FALSE)

  if(missing(fl)) {
    # we should rather parse the formula tree
    # find the terms involving G
    trm <- terms(formula,special='G')
    feidx <- attr(trm,'specials')$G-1
    festr <- paste(labels(trm)[feidx],collapse='+')

    if(festr == '') stop('No factors specified')
    # remove the G-terms from formula
    formula <- update(formula,paste('. ~ . -(',festr,')'))
    mf[['formula']] <- formula

    # then make a list of them, and find their names
    felist <- parse(text=paste('list(',gsub('+',',',festr,fixed=TRUE),')',sep=''))
    nm <- eval(felist,list(G=function(t) as.character(substitute(t))))
    # collapse them in case there's an interaction with a funny name
    nm <- lapply(nm,paste,collapse='.')
    
    # replace G with as.factor, eval with this, and the parent frame, or with data
    # allow interaction factors with '*'
    iact <- function(a,b) interaction(a,b,drop=TRUE)
    if(missing(data)) 
      fl <- eval(felist,list(G=as.factor,'*'=iact))
    else {
      G <- as.factor
      fl <- local({'*'<-iact;eval(felist,data,environment())})
    }

    names(fl) <- nm
  } else {
#    warning('The fl-argument is obsolete')
  }

  if(!is.list(fl)) stop('need at least one factor')
  if(!all(is.factor(unlist(fl,recursive=FALSE)))) {
    fl <- lapply(fl,as.factor)
  }
  if(is.null(names(fl))) names(fl) <- paste('fe',1:length(fl),sep='')

#  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

#  mf <- model.frame(formula,data=data,drop.unused.levels=TRUE)
  mt <- attr(mf,'terms')

  y <- model.response(mf,'numeric')

# try a sparse model matrix to save memory when removing intercept
# though, demeanlist must be full.  Ah, no, not much to save because
# it won't be sparse after centering
# we should rather let demeanlist remove the intercept, this
# will save memory by not copying.  But we need to remove it below in x %*% beta
# (or should we extend beta with a zero at the right place, it's only
#  a vector, eh, is it, do we not allow matrix lhs? No.)

  x <- model.matrix(mt,mf)
  rm(mf)
  icpt <- 0
  icpt <- which(attr(x,'assign') == 0)
  if(length(icpt) == 0) icpt <- 0
  ncov <- ncol(x) - (icpt > 0)
  if(ncov == 0) {
    # No covariates
    fr <- demeanlist(y,fl)
    z <- list(residuals=y,fe=fl,p=0,cfactor=compfactor(fl),full.residuals=fr,call=match.call())
    class(z) <- 'felm'
    return(z)
  }
  # here we need to demean things
  # in due time we should write demeanlist so that
  # it's sufficient with one call, so that y may be
  # parallelized together with the columns of x
  # (in a manner which avoids copying of data)
  dm <- demeanlist(list(y=y,x=x),fl,icpt)
#  cat(date(),'centering finished\n')
  yz <- dm[[1]]
  xz <- dm[[2]]
  rm(dm)
  badconv <- attr(xz,'badconv') + attr(yz,'badconv')
#  attributes(xz) <- attributes(x)
  dim(xz) <- c(nrow(x),ncov)
#  dimnames(xz) <- dimnames(x)
  attributes(yz) <- attributes(y)


# here we just do an lm.fit, however lm.fit is quite slow since
# it doesn't use blas (in particular it can't use e.g. threaded blas in acml)
# so we have rolled our own.

# we really don't return an 'lm' object or other similar stuff, so
# we should consider using more elementary operations which map to blas-3
# eg. solve(crossprod(xz),t(xz) %*% yz)
# Or, even invert by solve(crossprod(xz)) since we need
# the diagonal for standard errors.  We could use the cholesky inversion
# chol2inv(chol(crossprod(xz)))

  cp <- crossprod(xz)
  ch <- cholx(cp)

#  ch <- chol(cp)
#  beta <- drop(inv %*% (t(xz) %*% yz))
  # remove multicollinearities
  badvars <- attr(ch,'badvars')
  b <- crossprod(xz,yz)

  if(is.null(badvars)) {
    beta <- as.vector(backsolve(ch,backsolve(ch,b,transpose=TRUE)))
    inv <- chol2inv(ch)
  } else {
    beta <- rep(NaN,nrow(cp))
    beta[-badvars] <- backsolve(ch,backsolve(ch,b[-badvars],transpose=TRUE))
    inv <- matrix(NaN,nrow(cp),ncol(cp))
    inv[-badvars,-badvars] <- chol2inv(ch)
  }
  rm(b)
  if(icpt > 0) names(beta) <- colnames(x)[-icpt] else names(beta) <- colnames(x)
#  cat(date(),'projected system finished\n')
  z <- list(coefficients=beta,badconv=badconv)
  N <- nrow(xz)
  p <- ncol(xz) - length(badvars)

# how well would we fit with all the dummies?
# the residuals of the centered model equals the residuals
# of the full model, thus we may compute the fitted values
# resulting from the full model.
  zfit <- xz %*% ifelse(is.na(beta),0,beta)
  rm(xz)
  zresid <- yz - zfit
  rm(yz)
  z$fitted <- y - zresid
  z$full.residuals <- zresid
  # insert a zero at the intercept position
  if(icpt > 0) ibeta <- append(beta,0,after=icpt-1) else ibeta <- beta
  pred <- x %*% ifelse(is.na(ibeta),0,ibeta)
  z$residuals <- y - pred
#  z$xb <- pred
  rm(x)
  rm(y)
  gc()

  z$cfactor <- compfactor(fl)

  numrefs <- nlevels(z$cfactor) + max(length(fl)-2,0)
  numdum <- sum(unlist(lapply(fl,nlevels))) - numrefs
  z$numrefs <- numrefs
#  if(length(fl) <= 2) {
#    numdum <- sum(unlist(lapply(fl,nlevels))) - nlevels(z$cfactor)
#  } else {
#    numdum <- sum(unlist(lapply(fl,nlevels))) - length(fl) + 1
#  }
  z$df <- N - p - numdum
  vcvfactor <- sum(z$full.residuals**2)/z$df
  z$vcv <- inv * vcvfactor
  z$se <- sqrt(diag(z$vcv))
  z$sefactor <- sqrt(vcvfactor)
  z$tval <- z$coefficients/z$se
  z$pval <- 2*pt(abs(z$tval),z$df,lower.tail=FALSE)
  z$terms <- mt
  z$fe <- fl
  z$N <- N
  z$p <- p + numdum
  z$xp <- p
  z$call <- match.call()
  class(z) <- 'felm'

  return(z)
}

# return a data-frame with the group fixed effects, including zeros for references
getfe <- function(obj,references=NULL,se=FALSE,method='kaczmarz',ef='ref',bN=100) {

  if(method == 'kaczmarz') {
    if(!is.null(references))
       warning('use estimable function (ef) instead of references in the Kaczmarz method')
    if(is.null(ef)) ef <- 'ln'
    if(!is.character(ef) && !is.function(ef))
      stop('ef must be a function when using the Kaczmarz method')
    return(getfe.kaczmarz(obj,se,ef=ef,bN=bN))
  }
  if(method != 'cholesky') stop('method must be either kaczmarz or cholesky')
  attr(se,'sefactor') <- obj$sefactor
  attr(obj$fe,'references') <- references
  R <- obj$residuals
  # then the remaining.  This is usually sufficient.
  # we could also partition differently, just do the 'comp' adjustment accordingly
  # components
  ddlist <- makeddlist(obj$fe)
  gc()
  orignm <- attr(ddlist,'nm')
  comp <- 1
  res <- data.frame()
  for(dd in ddlist) {
#  res <- foreach(dd=ddlist,.combine=rbind,.init=data.frame()) %dopar% {
    dummies <- attr(dd,'dummies')
    keep <- attr(dd,'keep')
    comp <- attr(dd,'comp')
#    cat(date(),'comp dd',comp,'size',length(keep),'\n')
    Rhs <- as.vector(dummies %*% R[keep])
    names(Rhs) <- colnames(dd)
    alpha <- findfe(dd,Rhs,se)
    alpha[,'comp'] <- comp
    res <- rbind(res,alpha)
#    alpha
  }
  res <- res[orignm,]
  res[,'comp'] <- factor(res[,'comp'])    
# now, add factors telling which fe-group we're in
# the rownames are of the form <fe>.<idx>
  fefact <- strsplit(rownames(res),'.',fixed=TRUE)
  res[,'fe'] <- factor(unlist(lapply(fefact,function(l) l[[1]])))
  res[,'idx'] <- factor(unlist(lapply(fefact,function(l) paste(l[-1],collapse='.'))))
  return(res)
}

# Do a pivoted cholesky to detect multi-collinearities
cholx <- function(mat,eps=1e-4) {
  N <- dim(mat)[1]
  nm <- colnames(mat)
  o <- options(warn=-1)
  ch <- try(chol(mat),silent=TRUE)
  options(o)
  if(!inherits(ch,'try-error') && all(diag(ch) > eps*sqrt(diag(mat)))) return(ch)

  ch <- try(chol(mat,pivot=TRUE))
  pivot <- attr(ch,'pivot')

  rank <- min(which(c(diag(ch),-Inf) <= eps*c(sqrt(diag(mat)[pivot]),1)))-1
  okcol=1:rank
  badvars <- sort(pivot[-okcol])
  ok <- sort(pivot[okcol])
  return(structure(chol(mat[ok,ok]),badvars=badvars))
}
