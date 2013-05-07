# Author: Simen Gaure
# Copyright: 2011, Simen Gaure
# Licence: Artistic 2.0

.onLoad <- function(libname,pkgname) {
  options(lfe.usecg=FALSE)
  if(is.null(getOption('lfe.eps')))
    options(lfe.eps=1e-8)
  if(is.null(getOption('lfe.pint')))
    options(lfe.pint=300)
  if(is.null(getOption('lfe.threads'))) {
    cr <- as.integer(Sys.getenv('LFE_THREADS'))
    if(is.na(cr)) cr <- as.integer(Sys.getenv('OMP_NUM_THREADS'))
    if(is.na(cr)) cr <- numcores()
    if(is.na(cr)) {
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
  res <- .Call(C_demeanlist,
     mtx,
     as.list(fl),
     as.integer(icpt),               
     as.double(eps),
     as.integer(threads),
     as.integer(progress))

  if(!islist) res <- res[[1]]
  names(res) <- names(mtx)
  return(res)
}

compfactor <- function(fl) {
  if(length(fl) == 1) return(factor(rep(1,length(fl[[1]]))))
  factor(.Call(C_conncomp,fl[1:2]))
}


# a function for group fixed effect lm, i.e. first complete centering
# on category variables, then ols

flfromformula <- function(formula,data) {

  trm <- terms(formula,special='G')
  feidx <- attr(trm,'specials')$G-1
  festr <- paste(labels(trm)[feidx],collapse='+')

#  print(formula)
  if(festr == '') return(list(formula=formula,fl=NULL))
  # remove the G-terms from formula
  formula <- update(formula,paste('. ~ . -(',festr,')'))

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
  gc()
  names(fl) <- nm
  fl <- lapply(fl,as.factor)
  if(is.null(names(fl))) names(fl) <- paste('fe',1:length(fl),sep='')
  return(list(formula=formula,fl=fl))
}

# compute weighted cross product of M, i.e. where each row of M
# is scaled by the corresponding weight in w
# we're a bit lazy here, using a copy of the 
# data matrix M (note that we rely on w being recycled for each column)
wcrossprod <- function(M,w) {
  return(crossprod(M*as.vector(w)))
}

ccrossprod <- function(M,w,f) {
  # This is similar to the wcrossprod.
  # for level j, let u_j = sum(i in level j) w_i*x_i
  # then, compute the sum over the levels sum u_j' * u_j

  Mw <- M*as.vector(w)
  cM <- rowsum(Mw,f)
  return(crossprod(cM))
}

# ivresid is optional, used in 2. stage of 2sls to pass
# the difference between the original endogenous variable and the prediction
# for the purpose of computing sum of square residuals
doprojols <- function(psys, ivresid=NULL, exactDOF=FALSE) {
  if(is.null(psys$xz)) {
    # No covariates
    z <- list(residuals=psys$y,fe=psys$fl,p=0,cfactor=compfactor(psys$fl),full.residuals=psys$yz,call=match.call())
    class(z) <- 'felm'
    return(z)

  }
  yz <- psys$yz
  xz <- psys$xz
  y <- psys$y
  x <- psys$x
  fl <- psys$fl
  icpt <- psys$icpt
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
  z <- list(coefficients=beta,badconv=psys$badconv)
  z$N <- nrow(xz)
  z$p <- ncol(xz) - length(badvars)
  z$inv <- inv

# how well would we fit with all the dummies?
# the residuals of the centered model equals the residuals
# of the full model, thus we may compute the fitted values
# resulting from the full model.

# for the 2. step in the 2sls, we should replace 
# the instrumental variable with the real ones (in ivresid)
# when predicting, but only for the purpose of computing
# residuals.

  nabeta <- ifelse(is.na(beta),0,beta)
  zfit <- xz %*% nabeta

  zresid <- yz - zfit
  z$response <- y
  z$fitted <- y - zresid
  z$full.residuals <- zresid
  # insert a zero at the intercept position
  if(icpt > 0) ibeta <- append(beta,0,after=icpt-1) else ibeta <- beta

  pred <- x %*% ifelse(is.na(ibeta),0,ibeta)
  z$residuals <- y - pred

  z$ivresid <- rep(0,length(z$residuals))
  for(ivnam in names(ivresid)) {
#  if(!is.null(ivresid)) {
#    ivnam <- names(ivresid)[[1]]
    # fitted difference
    # predict the residuals from the first step
    ivbeta <- nabeta[[ivnam]]
    ivrp <- ivresid[[ivnam]] * ivbeta
    # subtract these to get the fit with the original variable
    z$full.residuals <- z$full.residuals - ivrp
    z$residuals <- z$residuals - ivrp
    z$ivresid <- z$ivresid + ivrp
  }

  z$terms <- psys$terms

  z$cfactor <- compfactor(fl)

  if(exactDOF && length(fl) > 2) {
    numrefs <- rankDefic(fl)
  } else {
    numrefs <- nlevels(z$cfactor) + max(length(fl)-2,0)
  }
  numdum <- sum(unlist(lapply(fl,nlevels))) - numrefs
  z$numrefs <- numrefs
  z$df <- z$N - z$p - numdum
  z$exactDOF <- exactDOF
  res <- z$full.residuals
  vcvfactor <- sum(res**2)/z$df
  z$vcv <- z$inv * vcvfactor
  
  # We should make the robust covariance matrix too.
  # it's inv * sum (X_i' u_i u_i' X_i) * inv
  # where u_i are the (full) residuals (Wooldridge, 10.5.4 (10.59))
  # i.e. inv * sum(u_i^2 X_i' X_i) * inv
  # for large datasets the sum is probably best computed by a series of scaled
  # rank k updates, i.e. the dsyrk blas routine, we make an R-version of it.
  # need to check this computation, the SE's are slightly numerically different from Stata's.
  # it seems stata does not do the small-sample adjustment
  dfadj <- (z$N-1)/z$df
  z$robustvcv <- dfadj * inv %*% wcrossprod(xz,res) %*% inv
  cluster <- psys$cluster

  # then the clustered covariance matrix, these do agree with stata, except that
  # stata does not do the degrees of freedom adjustment unless 'small' is specified
  if(!is.null(cluster)) {
    M <- nlevels(cluster)
    dfadj <- M/(M-1)*(z$N-1)/z$df
    z$clustervcv <- dfadj * inv %*% ccrossprod(xz,res,cluster) %*% inv
    z$cse <- sqrt(diag(z$clustervcv))
    z$ctval <- z$coefficients/z$cse
    z$cpval <- 2*pt(abs(z$ctval),z$df,lower.tail=FALSE)

  }

  z$se <- sqrt(diag(z$vcv))
  z$tval <- z$coefficients/z$se
  z$pval <- 2*pt(abs(z$tval),z$df,lower.tail=FALSE)

  z$rse <- sqrt(diag(z$robustvcv))
  z$rtval <- z$coefficients/z$rse
  z$rpval <- 2*pt(abs(z$rtval),z$df,lower.tail=FALSE)


  z$fe <- fl
# should we subtract 1 for an intercept?
# a similar adjustment is done in summary.felm when computing rdf
  z$p <- z$p + numdum - 1  
  z$xp <- z$p

  class(z) <- 'felm'
  return(z)
}

project <- function(mf,fl,data,clustervar=NULL) {

  if(!is.list(fl)) stop('need at least one factor')
  gc()


  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf,'terms')
  gc()
  if(!is.null(clustervar)) {
    if(!is.factor(clustervar)) {
      if(missing(data))
        clustervar <- as.factor(eval(as.name(clustervar),parent.frame()))
      else
        clustervar <- as.factor(data[,clustervar])
    }
  }
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
  gc()
  icpt <- 0
  icpt <- which(attr(x,'assign') == 0)
  if(length(icpt) == 0) icpt <- 0
  ncov <- ncol(x) - (icpt > 0)
  if(ncov == 0) {
    return(list(y=y,yz=demeanlist(list(y=y),fl=fl)[[1]],fl=fl))
  }

  # here we need to demean things
  dm <- demeanlist(list(y=y,x=x),fl,icpt)
  yz <- dm[[1]]
  xz <- dm[[2]]
  rm(dm)
  gc()

  badconv <- attr(xz,'badconv') + attr(yz,'badconv')
  dim(xz) <- c(nrow(x),ncov)
  if(icpt == 0)
    colnames(xz) <- colnames(x) 
  else 
    colnames(xz) <- colnames(x)[-icpt]
  attributes(yz) <- attributes(y)
  list(yz=yz,xz=xz,y=y,x=x,fl=fl,icpt=icpt,badconv=badconv,terms=mt,cluster=clustervar)  
}

felm <- function(formula, data, iv=NULL, clustervar=NULL, exactDOF=FALSE) {
  mf <- match.call(expand.dots = FALSE)


# If it's an iv-model, we should do as follows:
# take the 2. stage formula in 'formula', replace the
# rhs with the rhs of the iv-formula, remove this variable
# from the lhs and add the lhs of the iv-formula

  if(is.null(iv)) {
    # no iv, just do the thing
    sep <- flfromformula(formula,data)
    fl <- sep[['fl']]
    formula <- sep[['formula']]
    mf[['formula']] <- formula
    
    psys <- project(mf,fl,data,clustervar)
    z <- doprojols(psys,exactDOF=exactDOF)
    rm(psys)
    gc()
    z$call <- match.call()
    return(z)
  }

  # we must do the 1. step for each instrumented variable
  origform <- formula
  if(!is.list(iv)) iv <- list(iv)
  # collect the instrumented variables and remove them from origform
  baseform <- formula
  for(ivv in iv) {
    ivnam <- ivv[[2]]
  # create the new formula
    baseform <- update(baseform, substitute(. ~ . - Z,list(Z=ivnam)))
  }

  # then do the sequence of 1. steps
  # we put the instrumented variable on the lhs, and add in the formula for it on the rhs
  step2form <- origform
  ivarg <- list()
  vars <- NULL
  for(ivv in iv) {
     formula <- update(baseform, substitute(Z ~ . + V,list(Z=ivv[[2]],V=ivv[[3]])))
     sep <- flfromformula(formula,data)
     fl <- sep[['fl']]
     formula <- sep[['formula']]
     mf[['formula']] <- formula
     
     psys <- project(mf,fl,data,clustervar)
     z <- doprojols(psys)
     rm(psys)
     gc()
     # then we lift the fitted variable and create a new name
     ivz <- z
     evar <- as.character(ivv[[2]])
     new.var <- paste(evar,'(fit)',sep='')
     assign(new.var,ivz$fitted,envir=parent.frame())
     vars <- c(vars,new.var)
     # keep the residuals, they are need in doprojols below
     ivarg[[paste('`',new.var,'`',sep='')]] <- ivz$full.residuals
     # and add it to the equation
     step2form <- update(step2form,substitute(. ~ X - IV + FIT,
                                         list(X=step2form[[3]],IV=ivv[[2]],FIT=as.name(new.var))))

  }


  # now we have a formula in step2form with all the iv-variables
  # it's just to project it
  sep <- flfromformula(step2form,data)
  fl <- sep[['fl']]
  formula <- sep[['formula']]
  mf[['formula']] <- formula
  psys <- project(mf,fl,data,clustervar)
  rm(list=vars,envir=parent.frame())
  z <- doprojols(psys,ivresid=ivarg,exactDOF=exactDOF)
  rm(psys)
#  z$ivz <- ivz
  z$call <- match.call()
  return(z)
}

felm.old <- function(formula,fl,data) {
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
    gc()
    names(fl) <- nm
  } else {
#    warning('The fl-argument is obsolete')
  }

  if(!is.list(fl)) stop('need at least one factor')
  fl <- lapply(fl,as.factor)
  if(is.null(names(fl))) names(fl) <- paste('fe',1:length(fl),sep='')

  gc()
#  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf,'terms')
  gc()
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
  gc()
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

  dm <- demeanlist(list(y=y,x=x),fl,icpt)
  yz <- dm[[1]]
  xz <- dm[[2]]
  rm(dm)
  gc()

  badconv <- attr(xz,'badconv') + attr(yz,'badconv')
  dim(xz) <- c(nrow(x),ncov)
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
