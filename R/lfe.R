# Author: Simen Gaure
# Copyright: 2011, Simen Gaure
# Licence: Artistic 2.0

.onLoad <- function(libname,pkgname) {
  options(lfe.usecg=FALSE)
  if(is.null(getOption('lfe.eps')))
    options(lfe.eps=1e-8)
  if(is.null(getOption('lfe.pint')))
    options(lfe.pint=300)
  if(is.null(getOption('lfe.accel')))
    options(lfe.accel=1)
  if(is.null(getOption('lfe.bootmem')))
    options(lfe.bootmem=500)
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
		       progress=getOption('lfe.pint'),
                       accel=getOption('lfe.accel'),
                       randfact=TRUE) {

  if(length(fl) == 0) return(mtx)
  if(is.null(threads)) threads <- 1
  islist <- is.list(mtx)
  if(!islist) mtx <- list(mtx)
# randomize factor order, this may improve convergence
  if(randfact && length(fl) > 2) fl <- fl[order(runif(length(fl)))]

  res <- .Call(C_demeanlist,
     mtx,
     as.list(fl),
     as.integer(icpt),               
     as.double(eps),
     as.integer(threads),
     as.integer(progress),
     as.integer(accel))

  if(!islist) {
    res <- res[[1]]
    names(res) <- names(mtx)
  }
  res
}

edemeanlist <- function(...,fl,icpt=0,eps=getOption('lfe.eps'),
                       threads=getOption('lfe.threads'),
		       progress=getOption('lfe.pint'),
                       accel=getOption('lfe.accel'),
                       randfact=TRUE) {

  if(length(fl) == 0) return(list(...))
  if(is.null(threads)) threads <- 1
# randomize factor order, this may improve convergence
  if(randfact && length(fl) > 2) fl <- fl[order(runif(length(fl)))]

  .External(C_edemeanlist,
     ...,
     fl=as.list(fl),
     icpt=as.integer(icpt),               
     eps=as.double(eps),
     threads=as.integer(threads),
     progress=as.integer(progress),
     accel=as.integer(accel))
}
 
compfactor <- function(fl, WW=FALSE) {
  purefls <- sapply(fl,function(f) is.null(attr(f,'x')))
  fl <- fl[purefls]
  if(length(fl) <= 1) return(factor(rep(1,length(fl[[1]]))))
  if(WW && length(fl) > 2) {
    cf <- factor(.Call(C_wwcomp,fl))
  } else {
    cf <- factor(.Call(C_conncomp,fl[1:2]))
  }
  cf
}


# a function for group fixed effect lm, i.e. first complete centering
# on category variables, then ols

# parse our formula
flfromformula <- function(formula,data) {

  trm <- terms(formula,special='G')
  feidx <- attr(trm,'specials')$G+1
  va <- attr(trm,'variables')
  festr <- paste(sapply(feidx,function(i) deparse(va[[i]])),collapse='+')

  if(festr == '') return(list(formula=formula,fl=NULL))
  # remove the G-terms from formula
  formula <- update(formula,paste('. ~ . -(',festr,')'))

    # then make a list of them, and find their names
  felist <- parse(text=paste('list(',gsub('+',',',festr,fixed=TRUE),')',sep=''))
  nm <- eval(felist,list(G=function(arg) deparse(substitute(arg))))

#function(arg) {
#    deparse(substitute(arg))
#  }))
  # replace G with as.factor, eval with this, and the parent frame, or with data
  # allow interaction factors with '*'
  Gfunc <- function(f) if(is.null(attr(f,'xnam'))) factor(f) else f
  Ginfunc <- function(x,f) {
    if(is.factor(x))
      structure(interaction(factor(f),factor(x),drop=TRUE),xnam=deparse(substitute(x)),fnam=deparse(substitute(f)))
    else
      structure(factor(f),x=x,xnam=deparse(substitute(x)), fnam=deparse(substitute(f)))
  }

  if(missing(data)) 
    fl <- eval(felist,list(G=Gfunc, ':'=Ginfunc))
  else {
    G <- Gfunc
    fl <- local({':'<-Ginfunc; eval(felist,data,environment())})
  }

  names(fl) <- nm
#  cat('nm:\n')
#  print(unlist(nm))
#  cat('new nm:\n')
#  print(unlist(sapply(nm,function(n) if(is.null(attr(fl[[n]],'x'))) n else attr(fl[[n]],'fnam'))))
#  names(fl) <- sapply(nm,function(n) if(is.null(attr(fl[[n]],'x'))) n else attr(fl[[n]],'fnam'))


#  fl <- lapply(fl,as.factor)
  if(is.null(names(fl))) names(fl) <- paste('fe',1:length(fl),sep='')
  return(list(formula=formula,fl=fl))
}

# parse 
# use 2-part Formulas without G() syntax, like
# y ~ x1 + x2 | f1+f2
# or 3-part or more Formulas with iv-specification like
# y ~ x1 + x2 | f1+f2 | (q+w ~ x3+x4) | (z+v ~x5+x6)

## fl3fromformula <- function(form, data) {
##   f <- as.Formula(form)
##   opart <- formula(f,rhs=1)
##   gpart <- formula(f,lhs=0,rhs=2)
##   ivpart <- NULL
##   if(length(f)[[2]] > 2) {
##     ividx <- 3:(length(f)[[2]])
##     ivparts <- as.Formula(formula(f,lhs=0,rhs=ividx))
##   # Now, make a list of the iv-formulas where we split the lhs in each
##   # to obtain q ~ x3+x4, w ~x3+x4, z~x5+x6, v ~x5+x6
##     numivp <- length(ivparts)[2]
##     ivpart <- lapply(1:numivp,
##            function(i) {
##              ivf <- as.Formula(formula(ivparts,rhs=i))
##              # here we have something like q+w ~x3+x4
##              # loop through the lhs variables and make formulas
##              ff <- as.Formula(ivf[[2]][[2]])
##              lh <- formula(ff,rhs=0,drop=FALSE)
##              lapply(all.vars(lh), function(v) {
##                # we update the lhs of ff with v
##                update(ff,as.formula(substitute(L ~ .,list(L=as.name(v)))))
##              })
##            })
##   }
##   list(opart=opart,gpart=gpart,ivpart=unlist(ivpart,FALSE))
## }



# compute weighted cross product of M, i.e. where each row of M
# is scaled by the corresponding weight in w
# we're a bit lazy here, using a copy of the 
# data matrix M (note that we rely on w being recycled for each column)
wcrossprod <- function(M,w) {
  crossprod(M*as.vector(w))
}

ccrossprod <- function(M,w,f) {
  # This is similar to the wcrossprod.
  # for level j, let u_j = sum(i in level j) w_i*x_i
  # then, compute the sum over the levels sum u_j' * u_j

  Mw <- M*as.vector(w)
  cM <- rowsum(Mw,f)
  rm(Mw); gc()
  crossprod(cM)
}

# ivresid is optional, used in 2. stage of 2sls to pass
# the difference between the original endogenous variable and the prediction
# for the purpose of computing sum of square residuals
doprojols <- function(psys, ivresid=NULL, exactDOF=FALSE) {
  if(is.null(psys$yxz$x)) {
    # No covariates
    z <- list(r.residuals=psys$y,fe=psys$fl,p=0,cfactor=compfactor(psys$fl),
              na.action=psys$na.action, contrasts=psys$contrasts,
              fitted.values=psys$y - psys$yxz$y,
              df=length(psys$y)-(totalpvar(psys$fl)-nrefs(psys$fl, psys$cfactor, exactDOF)),
              residuals=psys$yxz$y,call=match.call())
    class(z) <- 'felm'
    return(z)

  }
  yz <- psys$yxz$y
  xz <- psys$yxz$x
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
  b <- crossprod(xz,yz)
  if(TRUE) {
  ch <- cholx(cp)
#  ch <- chol(cp)
#  beta <- drop(inv %*% (t(xz) %*% yz))
  # remove multicollinearities
  badvars <- attr(ch,'badvars')

  if(is.null(badvars)) {
    beta <- as.vector(backsolve(ch,backsolve(ch,b,transpose=TRUE)))
    inv <- chol2inv(ch)
  } else {
    beta <- rep(NaN,nrow(cp))
    beta[-badvars] <- backsolve(ch,backsolve(ch,b[-badvars],transpose=TRUE))
    inv <- matrix(NA,nrow(cp),ncol(cp))
    inv[-badvars,-badvars] <- chol2inv(ch)
  }
  rm(ch)
  } else {
    # try svd pseudo-inverse instead
    inv <- pinv(cp)
    badvars <- attr(inv,'badvars')
    beta <- as.vector(inv %*% b)
    inv[badvars,] <- NA
    inv[,badvars] <- NA
  }
  rm(b,cp)
  gc()

  if(length(fl) > 0 && icpt > 0) 
    names(beta) <- colnames(x)[-icpt] else names(beta) <- colnames(x)

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
# the instrumental variable with the real ones (the difference is in ivresid)
# when predicting, but only for the purpose of computing
# residuals.

  nabeta <- ifelse(is.na(beta),0,beta)
  zfit <- xz %*% nabeta
  zresid <- yz - zfit

  z$response <- y
  z$fitted.values <- y - zresid
  z$residuals <- zresid
  z$contrasts <- psys$contrasts
  if(length(fl) > 0) {
  # insert a zero at the intercept position
    if(icpt > 0) ibeta <- append(beta,0,after=icpt-1) else ibeta <- beta
    pred <- x %*% ifelse(is.na(ibeta),0,ibeta)
    z$r.residuals <- y - pred
  } else {
    z$r.residuals <- zresid
  }
  rm(x)

  gc()
  
  if(!is.null(ivresid)) {
    z$ivresid <- rep(0,length(z$r.residuals))
    for(ivnam in names(ivresid)) {
    # fitted difference
    # predict the residuals from the first step
      ivbeta <- nabeta[[ivnam]]
      ivrp <- ivresid[[ivnam]] * ivbeta
    # use this to subtract from residuals to get the residuals with the original variable
      z$ivresid <- z$ivresid + ivrp
    }
  }
  
  z$terms <- psys$terms

  z$cfactor <- compfactor(fl)

  totlev <- totalpvar(fl)

  if(is.numeric(exactDOF)) {
    z$df <- exactDOF
    numdum <- z$N - z$p - z$df
    z$numrefs <- totlev - numdum
  } else {
    numrefs <- nrefs(fl, z$cfactor, exactDOF)
    numdum <- totlev - numrefs
    z$numrefs <- numrefs
    z$df <- z$N - z$p - numdum
  }

  z$exactDOF <- exactDOF
  res <- z$residuals
  if(!is.null(ivresid)) res <- res - z$ivresid
  vcvfactor <- sum(res**2)/z$df
  z$vcv <- z$inv * vcvfactor

  setdimnames(z$vcv, list(names(beta),names(beta)))
#  dimnames(z$vcv) <- list(names(beta),names(beta))

  # We should make the robust covariance matrix too.
  # it's inv * sum (X_i' u_i u_i' X_i) * inv
  # where u_i are the (full) residuals (Wooldridge, 10.5.4 (10.59))
  # i.e. inv * sum(u_i^2 X_i' X_i) * inv
  # for large datasets the sum is probably best computed by a series of scaled
  # rank k updates, i.e. the dsyrk blas routine, we make an R-version of it.
  # need to check this computation, the SE's are slightly numerically different from Stata's.
  # it seems stata does not do the small-sample adjustment
  dfadj <- (z$N)/z$df
  # Now, here's an optimzation for very large xz. If we use the wcrossprod and ccrossprod
  # functions, we can't get rid of xz, we end up with a copy of it which blows away memory.
  # we need to scale xz with the residuals in xz, but we don't want to expand res to a full matrix,
  # and even get a copy in the result.
  # thus we modify it in place with a .Call. The scaled variant is also used in the cluster computation.
#  z$robustvcv <- dfadj * inv %*% wcrossprod(xz,res) %*% inv

  .Call(C_scalecols, xz, res)
#  z$robustvcv <- dfadj * inv %*% crossprod(xz) %*% inv
  z$robustvcv <- dfadj * crossprod(xz %*% ifelse(is.na(inv),0,inv))
  setdimnames(z$robustvcv, dimnames(z$vcv))
#  dimnames(z$robustvcv) <- dimnames(z$vcv)
  cluster <- psys$cluster
  gc()
  # then the clustered covariance matrix, these do agree with stata, except that
  # stata does not do the degrees of freedom adjustment unless 'small' is specified
  if(!is.null(cluster)) {
    M <- nlevels(cluster)
    dfadj <- M/(M-1)*(z$N-1)/z$df

##    z$clustervcv <- dfadj * inv %*% ccrossprod(xz,res,cluster) %*% inv
    cM <- rowsum(xz,cluster)
    rm(xz); gc()
#    z$clustervcv <- dfadj * inv %*% crossprod(cM) %*% inv
    cc <- crossprod(cM)
    rm(cM); gc()
    tmp <- cc %*% inv
    rm(cc); gc()
    z$clustervcv <- dfadj * inv %*% tmp
    rm(tmp,inv); gc()

    setdimnames(z$clustervcv, dimnames(z$vcv))
#    dimnames(z$clustervcv) <- dimnames(z$vcv)
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
  z$na.action <- psys$na.action
  class(z) <- 'felm'
  z
}

project <- function(mf,fl,data,contrasts,clustervar=NULL,pf=NULL) {

#  if(!is.list(fl)) stop('need at least one factor')

  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  subspec <- mf[['subset']]
  mf <- eval(mf, pf)
  mt <- attr(mf,'terms')
  naact <- attr(mf,'na.action')
  if(!is.null(naact))
    naclass <- attr(naact,'na.omit')

  if(!is.null(clustervar)) {
    if(!is.factor(clustervar)) {
      if(missing(data))
        clustervar <- factor(eval(as.name(clustervar),pf))
      else
        clustervar <- factor(data[,clustervar])
    }
  }
  # we need to change clustervar and factor list to reflect
  # subsetting and na.action. na.action is ok, it's set as an attribute in mf
  # but subset must be done manually. It's done before na handling
  if(!is.null(subspec)) {
    subs <- eval(subspec,pf)
    clustervar <- clustervar[subs]
    fl <- lapply(fl, function(fac) factor(fac[subs]))
  }
  if(!is.null(naact)) {
    clustervar <- clustervar[-naact]
    fl <- lapply(fl,function(fac) factor(fac[-naact]))
  }

  ret <- list(fl=fl, na.action=naact,terms=mt,cluster=clustervar)  
  rm(mt,clustervar,naact)


  ret$y <- model.response(mf,'numeric')

# try a sparse model matrix to save memory when removing intercept
# though, demeanlist must be full.  Ah, no, not much to save because
# it won't be sparse after centering
# we should rather let demeanlist remove the intercept, this
# will save memory by not copying.  But we need to remove it below in x %*% beta
# (or should we extend beta with a zero at the right place, it's only
#  a vector, eh, is it, do we not allow matrix lhs? No.)

# we make some effort to avoid copying the data matrix below
# this includes assigning to lists in steps, with gc() here and there.
# It's done for R 3.0.2. The copy semantics could be changed in later versions.


  ret$x <- model.matrix(ret$terms,mf,contrasts)
  ret$contrasts <- attr(ret$x,'contrasts')

  icpt <- 0
  icpt <- which(attr(ret$x,'assign') == 0)
  if(length(icpt) == 0) icpt <- 0
  ret$icpt <- icpt
  
  ncov <- ncol(ret$x) - (icpt > 0)
  if(ncov == 0) {
    ret$x <- NULL
    ret$yxz <- list(y=demeanlist(ret$y,fl))
    ret$yx <- list(y=ret$y)
    return(ret)
#    return(list(yx=list(y=y),yxz=list(y=demeanlist(y,fl=fl)[[1]]),fl=fl,contrasts=contr,na.action=naact))
  }

  # here we need to demean things
  # we take some care so that unexpected copies don't occur
  # hmm, the list() copies the stuff. How can we avoid a copy
  # and still enable parallelization over y and x in demeanlist? A vararg demeanlist?
  # I.e. an .External version?
#  yx <- list(y=ret$y, x=ret$x)
#  gc()
#  ret$yxz <- demeanlist(yx,fl,icpt)
#  rm(fl,yx); gc()
  ret$yxz <- edemeanlist(y=ret$y,x=ret$x,fl=fl,icpt=icpt)
  ret$badconv <- attr(ret$yxz$x,'badconv') + attr(ret$yxz$y,'badconv')
  # use our homebrewn setdimnames instead of colnames. colnames copies.
  if(length(fl) > 0) {
    if(icpt == 0)
      setdimnames(ret$yxz$x, list(NULL,colnames(ret$x)))
    else 
      setdimnames(ret$yxz$x, list(NULL,colnames(ret$x)[-icpt]))
  }
  ret
}

# this call changes or adds the dimnames of obj without duplicating it
# it should only be used on objects with a single reference
# Our use of it is safe, and we don't export it.
setdimnames <- function(obj, nm) {
  .Call(C_setdimnames,obj,nm)
}
  

felm <- function(formula, data, iv=NULL, clustervar=NULL, exactDOF=FALSE, subset, na.action, contrasts=NULL) {

  mf <- match.call(expand.dots = FALSE)

  pf <- parent.frame()
  if(is.null(iv)) {
    # no iv, just do the thing
    sep <- flfromformula(formula,data)
    fl <- sep[['fl']]
    formula <- sep[['formula']]
    mf[['formula']] <- formula
    
    psys <- project(mf,fl,data,contrasts,clustervar,pf)
    gc()

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
  step1 <- list()
  mf[['iv']] <- NULL
  for(ivv in iv) {
    # Now, make the full instrumental formula, i.e. with the rhs expanded with the
    # instruments, and the lhs equal to the instrumented variable
     fformula <- update(baseform, substitute(Z ~ . + V,list(Z=ivv[[2]],V=ivv[[3]])))
     sep <- flfromformula(fformula,data)
     fl <- sep[['fl']]
     mf[['formula']] <- sep[['formula']]
     psys <- project(mf,fl,data,contrasts,clustervar,pf)
     gc()
     z <- doprojols(psys, exactDOF=exactDOF)
     mf[['formula']] <- fformula
     z$call <- mf
     rm(psys)
     gc()

     # now, we need an ftest between the first step with and without the instruments(null-model)
     # We need the residuals with and without the
     # instruments. We have them with the instruments, but must do another estimation
     # for the null-model

     mfnull <- mf
     nullform <- update(baseform, substitute(Z ~ ., list(Z=ivv[[2]])))
     sepnull <- flfromformula(nullform, data)
     mfnull[['formula']] <- sepnull[['formula']]
     znull <- doprojols(project(mfnull, sepnull[['fl']], data, contrasts, clustervar, pf),
                        exactDOF=exactDOF)
     

     z$iv1fstat <- ftest(z,znull)
     z$rob.iv1fstat <- ftest(z,znull,vcov=z$robustvcv)
     if(!is.null(clustervar))
       z$clu.iv1fstat <- ftest(z,znull,vcov=z$clustervcv)
     step1 <- c(step1,list(z))

     # then we lift the fitted variable and create a new name
     ivz <- z
     evar <- as.character(ivv[[2]])
     new.var <- paste(evar,'(fit)',sep='')
     assign(new.var,ivz$fitted.values,envir=parent.frame())
     vars <- c(vars,new.var)
     # keep the residuals, they are needed to reconstruct the residuals for the
     # original variables in the 2. stage
     ivarg[[paste('`',new.var,'`',sep='')]] <- ivz$residuals
     # and add it to the equation
     step2form <- update(step2form,substitute(. ~ X - IV + FIT,
                                         list(X=step2form[[3]],IV=ivv[[2]],FIT=as.name(new.var))))
     
  }
  names(step1) <- names(iv)

  # now we have a formula in step2form with all the iv-variables
  # it's just to project it
  sep <- flfromformula(step2form,data)
  fl <- sep[['fl']]
  formula <- sep[['formula']]
  mf[['formula']] <- formula
  psys <- project(mf,fl,data,contrasts,clustervar,pf)
  rm(list=vars,envir=parent.frame())
  z <- doprojols(psys,ivresid=ivarg,exactDOF=exactDOF)
  z$step1 <- step1
  rm(psys)
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

#  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
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
    z <- list(r.residuals=y,fe=fl,p=0,cfactor=compfactor(fl),residuals=fr,call=match.call())
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
  z$fitted.values <- y - zresid
  z$residuals <- zresid
  # insert a zero at the intercept position
  if(length(fl) > 0) {
    if(icpt > 0) ibeta <- append(beta,0,after=icpt-1) else ibeta <- beta
    pred <- x %*% ifelse(is.na(ibeta),0,ibeta)
    z$r.residuals <- y - pred
  } else {
    z$r.residuals <- zresid
  }
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
  vcvfactor <- sum(z$residuals**2)/z$df
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
getfe <- function(obj,references=NULL,se=FALSE,method='kaczmarz',ef='ref',bN=100, robust=FALSE, cluster=NULL) {

  if(length(obj$fe) == 0) return(NULL)
  if(method == 'kaczmarz') {
    if(!is.null(references))
       warning('use estimable function (ef) instead of references in the Kaczmarz method')
    if(is.null(ef)) ef <- 'ln'
    if(!is.character(ef) && !is.function(ef))
      stop('ef must be a function when using the Kaczmarz method')
    return(getfe.kaczmarz(obj,se,ef=ef,bN=bN, robust=robust, cluster=cluster))
  }
  if(method != 'cholesky') stop('method must be either kaczmarz or cholesky')
  attr(se,'sefactor') <- obj$sefactor
  attr(obj$fe,'references') <- references
  R <- obj$r.residuals
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

