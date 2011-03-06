# Author: Simen Gaure
# Copyright: 2011, Simen Gaure
# Licence: Artistic 2.0

.onLoad <- function(libname,pkgname) {
  registerDoSEQ()
  o <- options(warn=-1)
  require('doMC',quietly=TRUE,warn.conflicts=FALSE)
  if(!exists('registerDoMC')) {
    require('doSMP',quietly=TRUE,warn.conflicts=FALSE)
  }
  options(o)
  usePAR <- exists('registerDoMC') || exists('registerDoSMP')
  if(is.null(getOption('lfe.eps')))
    options(lfe.eps=1e-8)
  if(is.null(getOption('lfe.threads'))) {
    cr <- as.integer(Sys.getenv('LFE_THREADS'))
    if(is.na(cr)) cr <- as.integer(Sys.getenv('OMP_NUM_THREADS'))
    if(is.na(cr)) cr <- 1
    options(lfe.threads=cr)
    if(cr == 1) cat('LFE using 1 thread, set options(lfe.threads=n) to use more\n')
  }

  if(!usePAR)
    message('on multicore machines, bootstrapping will run faster if package "doMC" or "doSMP" is installed')

}

#   Some things in this file is done in a weird way.
#   In some cases there are efficiency reasons for this, e.g. because
#   the "standard" way of doing things may result in a copy which is costly
#   when the problem is *large*.
#   In other cases it may simply be due to the author's unfamiliarity with how
#   things should be done in R


demeanlist <- function(mtx,fl,icpt=0,eps=getOption('lfe.eps'),
                       threads=getOption('lfe.threads')) {
  if(is.null(threads)) threads <- 1
  islist <- is.list(mtx)
  if(!islist) mtx <- list(mtx)
  res <- .Call('demeanlist',
     mtx,
     as.list(fl),
     as.integer(icpt),               
     as.double(eps),
     as.integer(threads),PACKAGE='lfe')

  if(!islist) res <- res[[1]]
  names(res) <- names(mtx)
  return(res)
}

compfactor <- function(fl) {
  if(length(fl) == 1) return(factor(rep(1,length(fl[[1]]))))
  factor(.Call('conncomp',fl[1:2]))
}

kaczmarz <- function(fl,R,eps=getOption('lfe.eps')) {
  v <- .Call('kaczmarz',fl,as.vector(R),eps)
  names(v) <- unlist(lapply(names(fl),function(n) paste(n,levels(fl[[n]]),sep='.')))
  v
}

# a function for removing the intercept from a model matrix
# lifted from package 'pls' and modified.  Not used at the moment
delete.intercept <- function (mm) 
{
    saveattr <- attributes(mm)
    intercept <- which(saveattr$assign == 0)
    saveattr$assign <- saveattr$assign[-intercept]
    if (!length(intercept)) 
        return(mm)
    mm <- mm[, -intercept, drop = FALSE]
    for(a in names(saveattr)) {
      
      if(inherits(mm,'Matrix')) {
        if(a %in% c('x','Dim','Dimnames','factors','class')) next      
      } else {
        if(a %in% c('dim','dimnames')) next
      }
      attr(mm,a) <- saveattr[[a]]
    }
    mm
}



findfe <- function(dd,Rhs,se=FALSE) {
  # find references
  refnames <- attr(dd,'refnames')
  nm <- c(colnames(dd),refnames)
  refcnt <- attr(dd,'refcnt')
  # add in the reference in dg

  dg <- c(diag(dd),refcnt)
  ok <- 1:ncol(dd)

  if(se) sev <- double(length(nm))
  alphacoef <- double(length(nm))

  # the super-nodal algorithm
  # is default and far better, but it consumes more memory

  trysolve <- try(solve(dd,Rhs))
  if(inherits(trysolve,'try-error')) {
    if(grepl('problem too large',geterrmessage())) {
      message(paste('Never mind, trying *slower* non-supernodal algorithm, nnz=',nnzero(dd)))
      message(paste(date(),'This may be an opportunity for a nice cup of tea. Or two.'))
      gc()
      ch <- Cholesky(dd,super=FALSE,perm=TRUE)
      trysolve <- solve(ch,Rhs)
      rm(ch)
      gc()
    } else {
      stop(geterrmessage())
    }
  } 


  alphacoef[ok] <- as.vector(trysolve)
  if(se) {
    # is there a faster way to find the diagonal of the inverse?
    sev[ok] <- sqrt(diag(solve(dd)))*attr(se,'sefactor')
    alpha <- data.frame(effect=alphacoef,se=sev,obs=dg)
  } else {
    alpha <- data.frame(effect=alphacoef,obs=dg)
  }
  rownames(alpha) <- nm
  alpha
}


makedummies <- function(factors) {
  nm <- c()
  dummies <- Matrix(0,0,length(factors[[1]]))
  for(i in 1:length(factors)) {
    f <- factors[[i]]
    dummies <- rBind(dummies,as(f,'sparseMatrix'))
    nm <- c(nm,paste(names(factors)[[i]],levels(f),sep='.'))
  }
  rownames(dummies) <- nm
  dummies
}

makedd.full <- function(factors) {
  dm <- makedummies(factors)
  nm <- rownames(dm)
  dd <- tcrossprod(dm)
  rownames(dd) <- colnames(dd) <- nm
  attr(dd,'dummies') <- dm
  attr(dd,'nm') <- nm
  dd
}

makeddlist <- function(factors) {
  if(length(factors) > 2) {
    if(is.null(attr(factors,'references'))) {
      # find references by fiddling with Cholesky
      message('*** More than two groups, finding refs by Cholesky pivots, interpret at own risk')
      # first the full matrix, find small pivots
      dd <- makedd.full(factors)
      orignm <- attr(dd,'nm')

      # add small amount to diagonal
      eps <- sqrt(.Machine$double.eps)
      Ch <- try(Cholesky(dd,super=TRUE,perm=TRUE,Imult=eps))
      if(inherits(Ch,'try-error') && grepl('problem too large',geterrmessage())) {
        Ch <- Cholesky(dd,super=FALSE,perm=TRUE,Imult=eps)
      }
      # strangely enough, coercing to sparseMatrix doesn't take care of
      # the permutation, we apply it manually.  Let's hope it's never fixed.
      rm(dd); gc()
      pivot <- Ch@perm
      ch <- as(Ch,'sparseMatrix')
      rm(Ch); gc()
      dg <- diag(ch)[order(pivot)]**2
      rm(ch); gc()
      refs <- (dg < eps**(1/3))
      refnames <- orignm[refs]
      message(paste('***',length(refnames),'references found'))
    } else {
      refnames <- attr(factors,'references')
      orignm <- unlist(lapply(names(factors),
                              function(n) paste(n,levels(factors[[n]]),sep='.')))
    }
    # there may be references in more than one factor
    # remove all of them
    # create factor list with named levels
    nf <- lapply(names(factors),function(n) {
      f <- factors[[n]]
      levels(f) <- paste(n,levels(f),sep='.')
      f
    })
    # remove reference levels, and remove the prefix
    # find the levels
    lev <- lapply(nf,function(f) which(levels(f) %in% refnames))
    nnf <- mapply(function(f,l) factor(f,exclude=levels(f)[l]),factors,lev,SIMPLIFY=FALSE)
    dd <- makedd.full(nnf)
    attr(dd,'keep') <- 1:length(nnf[[1]])
    attr(dd,'refnames') <- refnames
#    attr(dd,'refcnt') <- rep(1,length(refnames))
# find the number of occurences
    cntlst <- unlist(lapply(refnames,function(n) lapply(nf,function(f) sum(f == n))))
    attr(dd,'refcnt') <- cntlst[cntlst > 0]
    attr(dd,'comp') <- 1
    res <- list(dd)
    attr(res,'nm') <- orignm
  } else {
    # 2 or fewer factors, find references by component
    cf <- compfactor(factors)
    nml <- lapply(factors,function(f) levels(f))
    nm <- unlist(lapply(names(nml),function(n) paste(n,nml[[n]],sep='.')))
    res <- list()
    li <- 1
# this loop suffers from too much copying and stuff
# when there are many components (e.g. like 10000)
    remfact <- factors
    fullidx <- 1:length(factors[[1]])
    for(l in levels(cf)) {
      # find those in this level
      keep <- which(cf == l)
#      cat(date(),'comp fact',li,'size',length(keep),'\n')
      fcomp <- lapply(remfact,function(f) factor(f[keep]))
      remfact <- lapply(remfact,function(f) factor(f[-keep]))
      cf <- factor(cf[-keep])
      # then the reference level
      maxrefs <- lapply(fcomp,function(f) {tf <- table(f); m <- which.max(tf); tf[m]})
      # in which factor
      rfac <- which.max(unlist(maxrefs))
      # which level
      reflevel <- names(maxrefs[[rfac]])
      # drop that level from the factor
      fcomp[[rfac]] <- factor(fcomp[[rfac]],exclude=reflevel)
      refname <- paste(names(remfact)[[rfac]],reflevel,sep='.')
      # remove those without levels
      len <- unlist(lapply(fcomp,nlevels))
      fcomp <- fcomp[len > 0]
      dd <- makedd.full(fcomp)
      # the keep attribute should be relative to
      # the full factor, not to remfact
      attr(dd,'keep') <- fullidx[keep]
      fullidx <- fullidx[-keep]
      attr(dd,'refnames') <- refname
      attr(dd,'refcnt') <- max(unlist(maxrefs))
      attr(dd,'comp') <- li
      res[[li]] <- dd
      li <- li+1
#      res <- c(res, list(dd))
    }
    attr(res,'nm') <- nm
  }
  res
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
  z$residuals <- y - x %*% ifelse(is.na(ibeta),0,ibeta)
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

# A common estimable function on the fe-coefficients
# return an estimable function, the matrix which
# pins a reference in each component, the one with the
# most observations
# if there are more than two factors, assume they don't
# ruin identification beyond one extra reference for each such factor
defaultef <- function(est) {

# the names of the dummies, e.g. id.4 firm.23
  nm <- unlist(lapply(names(est$fe),function(n) paste(n,levels(est$fe[[n]]),sep='.')))

# how many obervations for each level
  obs <- unlist(lapply(est$fe,table))  

  if(length(est$fe) == 2) {
    comp <- unlist(lapply(est$fe, function(f) est$cfactor[match(levels(f),f)]))
  } else if(length(est$fe) > 2) {
    # we should formally assign unique component numbers for factors beyond the second
    comp <- factor(unlist(lapply(est$fe[1:2], function(f) est$cfactor[match(levels(f),f)])))
    exlvls <- (nlevels(comp)+1):(nlevels(comp)+1 + length(est$fe)-3)
    comp <- c(comp,mapply(rep,exlvls,unlist(lapply(est$fe[3:length(est$fe)],nlevels))))
  } else {
    comp <- rep(1,length(obs))
  }

  refnames <- unlist(tapply(obs,comp,function(l) names(which.max(l))))
  refno <- match(refnames,nm)
  refsub <- refno[comp]
  # now v[refnames] will be the reference values
  # we should do a v <-  v - v[refsub]
  # but for the main components we should only do this for the
  # factor in which the references is
  # for the other factor we should add the reference
  # thus we need two versions of refsub, one with NA's in the
  # reference factor, one with NA's in the other, then we must
  # replace NA's with zero before subtracting
  # so which ones belong to which factor?
  # make a factor to decide
  # we're now only interested in the two first factors, those which define the
  # components
  fef <- factor(unlist(lapply(names(est$fe),function(n) rep(n,nlevels(est$fe[[n]])))))
  # then figure out in which factor the reference is
  rf <- sub('(^.*)\\.[^.]*$','\\1',refnames)
  # now, create a refsubs which is the ones to be subtracted
  # each refsub belonging to somthing else than the reference factor
  # should be NA'ed.
  if(length(est$fe) > 2) {
    extra <- (length(refno)-length(est$fe)+3):length(refno)
    sw <- c(names(est$fe)[c(2,1)],rep('.NA',length(est$fe)-2))
  } else {
    swap <- if(length(est$fe) == 2) c(2,1) else 1
    sw <- names(est$fe)[swap]
    extra <- integer(0)
  }
  names(sw) <- names(est$fe)
  otherf <- sw[rf]
  # which should we keep subtracting?  Those which are different from the
  # the other factor, hmm, fix this test, rf is not a singleton
  # force these components to NA, which components?
  nosub <- fef != rf[comp]
  refsubs <- refsub
  refsubs[nosub] <- NA
  # which should we add, those which are different from the reference factor

  noadd <- fef != otherf[comp]
  refsuba <- refsub
  refsuba[noadd] <- NA
  extrarefs <- refno[extra]

  # now, what if we want zero-means on the other-factor?

  # return a function doing the stuff
  # set up a minimal environment for it

  local(function(v) {
    esum <- sum(v[extrarefs])
    df <- v[refsubs]
    sub <- ifelse(is.na(df),0,df)
    df <- v[refsuba]
    add <- ifelse(is.na(df),0,df+esum)
    v - sub + add
  },list(extrarefs=extrarefs,refsubs=refsubs,refsuba=refsuba))
}

btrap <- function(alpha,obj,N=100,ef=NULL,eps=getOption('lfe.eps'),threads=getOption('lfe.threads')) {
  # bootstrap the stuff
  # bootstrapping is really to draw residuals over again, i.e. to change
  # the outcome.  Predictions of the estimated system are adjusted by
  # drawing from the residuals. We do have PY-PXbeta, we resample it,
  # call it nres, then let newR = (I-P)*nres, will this do the trick?
  # it ought to do, but it's the variation conditional on the least norm
  # solution, it's not directly relevant to effects.

  if(is.null(ef)) 
    ef <- attr(alpha,'ef')
  else {
    v <- ef(alpha[,'effect'])
    if(length(v) != nrow(alpha) || isTRUE(attr(ef,'noP'))) {
      alpha <- data.frame(effect=v)
      rownames(alpha) <- names(v)
    } else {
      alpha[,'effect'] <- v
    }
  }
  R <- obj$residuals-obj$full.residuals
  j <- 0 #avoid warning from check
  gc()
  if(threads > 1) {
    if(exists('registerDoMC')) {
      registerDoMC(threads)
    } else if(exists('registerDoSMP')) {
      workers <- startWorkers(threads)
      on.exit(stopWorkers(workers))
      registerDoSMP(workers)
      # load the dll for kaczmarz in the workers
      dll <- getNativeSymbolInfo('kaczmarz','lfe')[['package']][['name']]
      foreach(j=1:threads) %dopar% {
        library.dynam(dll,package='lfe')
      }
    } else {
      warning('No parallel backend available, try loading package "doMC" or "doSMP" to speed up bootstrapping')
    }
  } else {
    registerDoSEQ()
  }
  # get somewhat larger chunks
  thr <- 2*threads-1
  blk <- N %/% thr
  nsamp <- blk * thr
  rsamp <- matrix(0,length(R),thr)
  # do we really need Rthr?  Don't we get the statistically same result without?
  Rthr <- rep(R,thr)
  # we should really compute new residuals for new betas, now we assume
  # the beta's will be unchanged.  Is there some way we may incorporate
  # their variation without doing the estimation?
  sefact <- sqrt(var(obj$residuals)/var(obj$residuals-obj$full.residuals) * obj$df/obj$N)
  sm <- mean(obj$full.residuals) # ought to be zero
  smpdraw <- sefact*(as.vector(obj$full.residuals)-sm)#+sm

  vsum <- vector('double',nrow(alpha))
  vsq <-  vector('double',nrow(alpha))

  for(i in 1:blk) {
    rsamp[,] <- sample(smpdraw,length(R)*thr,replace=TRUE)
    newR <- rsamp - demeanlist(rsamp,obj$fe,eps=eps,threads=threads) + Rthr
    # we should have a parallel kaczmarz, we don't yet
    # so use foreach
    ke <- foreach::foreach(j=1:thr,.combine=cbind,.export=c('kaczmarz','ef')) %dopar% {
      ke <- kaczmarz(obj$fe,newR[,j],eps)
      # here we may contrast the coefficients, but how?
      # zero-mean, pinning to zero?
      # let user decide with an estimable function ef
      if(!is.null(ef)) ke <- ef(ke)
      ke
    }
    # in case thr==1:
    dim(ke) <- c(nrow(alpha),thr)
    vsum <- vsum + rowSums(ke)
    vsq <- vsq + rowSums(ke**2)
  }
  alpha[,'se'] <- sqrt(vsq/nsamp - (vsum/nsamp)**2)
  alpha
}

getfe.kaczmarz <- function(obj,se=FALSE,eps=getOption('lfe.eps'),ef=defaultef(obj),bN=100) {

  R <- obj$residuals-obj$full.residuals
  v <- kaczmarz(obj$fe,R,eps)
  nm <- names(v)
  if(!is.null(ef)) v <- ef(v)
  if(isTRUE(attr(ef,'noP')) || length(v) != length(nm)) {
    res <- data.frame(effect=v)
    rownames(res) <- names(v)
  } else {
    fegrp <- factor(unlist(lapply(names(obj$fe),function(n) rep(n,nlevels(obj$fe[[n]])))))
    feidx <- factor(unlist(lapply(names(obj$fe),function(n) levels(obj$fe[[n]]))))
    nm <- paste(fegrp,feidx,sep='.')
    obs <- unlist(lapply(obj$fe,table))
    if(length(obj$fe) == 2) {
      comp <- unlist(lapply(obj$fe, function(f) obj$cfactor[match(levels(f),f)]))
    } else if(length(obj$fe) > 2) {
      # we should formally assign unique component numbers for factors beyond the second
      comp <- factor(unlist(lapply(obj$fe[1:2], function(f) obj$cfactor[match(levels(f),f)])))
      exlvls <- (nlevels(comp)+1):(nlevels(comp)+1 + length(obj$fe)-3)
      comp <- c(comp,mapply(rep,exlvls,unlist(lapply(obj$fe[3:length(obj$fe)],nlevels))))
    } else {
      comp <- rep(1,length(obs))
    }
    res <- data.frame(effect=NA,se=NA,obs=obs,comp=factor(comp),fe=fegrp,idx=feidx)
    if(!se) res[,'se'] <- NULL

    res[,'effect'] <- v
    rownames(res) <- nm
  }

  attr(res,'ef') <- ef
  if(se) {
    res <- btrap(res,obj,bN,ef=ef,eps=eps)
  }
  res
}

# return a data-frame with the group fixed effects, including zeros for references
getfe <- function(obj,references=NULL,se=FALSE,method='kaczmarz',ef=defaultef(obj),bN=100) {

  if(method == 'kaczmarz') {
    if(!is.null(references))
       warning('use estimable function (ef) instead of references in the Kaczmarz method')
    if(!is.null(ef) && !is.function(ef))
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

# lfefile <- function(fn) {
#  base <- strsplit(getLoadedDLLs()[['lfe']][['path']],'/')[[1]]
#  scr <- paste(paste(base[1:(length(base)-2)],collapse='/'),'/exec/lfescript',sep='')
#  cmd <- paste('perl',scr,fn)
#  source(pipe(cmd,open='r'))
# }

# An idea I had to multi-fe's. Could we start out by picking two factors, find
# references there (assuming there are no more fe's.  We must at least have these
# references, don't we?).  Then replace these
# two factors with their interaction and do the same with this and the
# third factor? (of course not picking references in the interaction factor), and so on?
# I have no proof.  The first assertion is ok, but the reduction to fewer factors
# with the interaction is dubious.

# This algorithm is obviously not correct, the result (i.e. the number of refs)
# depends on the order of the factors.  I'll have to think more about it
# (Though, doing a regression on the covariates and the resulting effects,
# yields coefficients almost (but not exactly) equal to 1 for the fixed effects, as it should)
# It could be way more complicated, graph theory tends to be.

# this one is experimental, it tries to find the references
# when there are more than two factors.
# Only tested with individual,firm, year fixed effects.
# and individual,firm,shoe,year effects (some people change their shoes...)

## findrefs <- function(fl) {
##   stop("not working yet")
##   if(length(fl) == 1) {
##     return(levels(fl[[1]])[[which.max(rowSums(as(fl[[1]],'sparseMatrix')))]])
##   }
##   # order by decreasing number of levels
##   oo <- order(unlist(lapply(fl,nlevels)),decreasing=TRUE)
##   fl <- fl[oo]
##   refnames <- c()
##   first <- fl[[1]]
##   nextf <- 2
##   while(nextf <= length(fl)) {
##     f2 <- fl[[nextf]]
##     n2 <- names(fl)[[nextf]]
##     levels(f2) <- paste(n2,levels(f2),sep='.')
##     ft <- list(first=first)
##     ft[[n2]] <- f2
##     cf <- compfactor(ft)
##     # Now, within each component, find the f2 with most occurrences, this is a reference
##     # is there an easier way to do this?
##     ref <- unlist(tapply(f2,cf,
##                  function(fc) {
##                    fc <- factor(fc)
##                    levels(fc)[[which.max(rowSums(as(fc,'sparseMatrix')))]]
##                  }))
##     refnames <- c(refnames,ref)
##     if(nextf < length(fl)) {
##       # collapse the reference levels in f2 and drop unused levels
##       f2[f2 %in% ref] <- ref[[1]]
##       f2 <- factor(f2)
##       first <- factor(paste(first,f2,sep='.'))
##     }
##     nextf <- nextf+1
##   }
## # this may be too many references, how should we go about
## # removing unnecessary ones?
 
##   refnames
## }
