

kaczmarz <- function(fl,R,eps=getOption('lfe.eps'),init=NULL,threads=getOption('lfe.threads')) {
  if(getOption('lfe.usecg')) return(cgsol(fl,R,eps,init))
  if(is.null(threads)) threads <- 1
  islist <- is.list(R)
  if(!islist) R <- list(R)
  v <- .Call(C_kaczmarz,fl,R,eps,as.vector(init),as.integer(threads))
  if(!islist) {
    v <- drop(v[[1]])
  }
  v
}

cgsol <- function(fl,R,eps=getOption('lfe.eps'), init=NULL,threads=NULL) {
  if(is.list(R)) stop("cgsol can't handle list R")

  # create matrix
  mat <- as(fl[[1]],'sparseMatrix')
  for(f in fl[-1]) {
    mat <- rBind(mat,as(f,'sparseMatrix'))
  }
  mat <- t(mat)

# target function
  fn <- function(p) sum((mat %*% p - R)^2)
# its gradient
  gr <- function(p) {
      a <- 2*(mat%*%p-R)
      dim(a) <- c(1,nrow(mat))
      as.vector(a %*% mat)
  }
  if(is.null(init)) init <- rep(0,ncol(mat))
  cg <- Rcgmin(init,fn,gr,control=list(eps=eps))
  if(cg$convergence != 0) warning(cg$message)
  cg$par
}

getfe.kaczmarz <- function(obj,se=FALSE,eps=getOption('lfe.eps'),ef='ref',bN=100) {

  R <- obj$r.residuals-obj$residuals
  v <- kaczmarz(obj$fe,R,eps)
  nm <- names(v)
  if(is.character(ef)) {
    ef <- efactory(obj,opt=ef)
  }
  if(!isTRUE(attr(ef,'verified')) && !is.estimable(ef,obj$fe)) {
    warning('Supplied function seems non-estimable\n')
  }
  v <- ef(v,TRUE)
  res <- data.frame(effect=v)
  if(!is.null(attr(v,'extra'))) res <- cbind(res,attr(v,'extra'))
  rownames(res) <- names(v)
  attr(res,'ef') <- ef
  if(se) {
    res <- btrap(res,obj,bN,eps=eps)
  }
  res
}


# A common estimable function on the fe-coefficients
# return an estimable function, the matrix which
# pins a reference in each component, the one with the
# most observations
# if there are more than two factors, assume they don't
# ruin identification beyond one extra reference for each such factor

# Note: when I get some time, I'll implement a Weeks-Williams estimable
# function.  I.e. with a reference in each factor in each component
# from compfactor(...,WW=TRUE). This may be what many people want.

efactory <- function(obj, opt='ref', ...) {

# the names of the dummies, e.g. id.4 firm.23
  nm <- unlist(lapply(names(obj$fe),function(n) paste(n,levels(obj$fe[[n]]),sep='\003')))

# how many obervations for each level
  lobs <- lapply(obj$fe,table)
  obs <- unlist(lobs)  
  names(obs) <- unlist(lapply(names(lobs), function(n) paste(n,'\003',names(lobs[[n]]),sep='')))
  if(length(obj$fe) == 2) {
    # now, find the component of each parameter, i.e. each level.  We do this
    # by finding the first occurence of each level, i.e. match(levels(f),f)
    comp <- factor(unlist(lapply(obj$fe, function(f) obj$cfactor[match(levels(f),f)])))
    ncomp <- nlevels(comp)
  } else if(length(obj$fe) > 2) {
    # we should formally assign unique component numbers for factors beyond the second
    comp <- factor(unlist(lapply(obj$fe[1:2], function(f) obj$cfactor[match(levels(f),f)])))
    ncomp <- nlevels(comp)
    exlvls <- (nlevels(comp)+1):(nlevels(comp)+1 + length(obj$fe)-3)
    comp <- as.factor(c(comp,unlist(mapply(rep,exlvls,unlist(lapply(obj$fe[3:length(obj$fe)],nlevels))))))
  } else {
    comp <- factor(rep(1,length(obs)))
    ncomp <- 1
  }

  refnames <- unlist(tapply(obs,comp,function(l) names(which.max(l))))
  # now v[refnames] will be the reference values

  refno <- match(refnames,nm)
  refsub <- refno[comp]
  # refsub is a vector, in entry i it contains the reference for entry i
  # i.e. we should do a v <-  v - v[refsub]
  # subtract all references.
  # but for the main components we should only do this for the
  # factor in which the references is.
  # for the other factor we should add the reference
  # thus we need two versions of refsub, one with NA's in the
  # reference factor, one with NA's in the other, then we must
  # replace NA's with zero before subtracting
  # so which ones belong to which factor?
  # make a factor to decide

  fef <- factor(unlist(lapply(names(obj$fe),function(n) rep(n,nlevels(obj$fe[[n]])))))
  # level of the factor
  idx <- factor(unlist(lapply(obj$fe,function(f) levels(f))))
  # then figure out in which factor the reference is
  # make sure to allow '.' in factor names
  rf <- sub('(^.*)\003..*$','\\1',refnames)
  # now, create a refsubs which is the ones to be subtracted
  # each refsub belonging to somthing else than the reference factor
  # should be NA'ed.
  if(length(obj$fe) > 2) {
    extra <- (length(refno)-length(obj$fe)+3):length(refno)
    sw <- c(names(obj$fe)[c(2,1)],rep('.NA',length(obj$fe)-2))
  } else {
    swap <- if(length(obj$fe) == 2) c(2,1) else 1
    sw <- names(obj$fe)[swap]
    extra <- integer(0)
  }
  names(sw) <- names(obj$fe)
  otherf <- sw[rf]
  # which should we not subtract?  
  # Those which are 
  nosub <- fef != rf[comp]
  refsubs <- refsub
  refsubs[nosub] <- NA

  # which should we add, those which are different from the reference factor
  noadd <- fef != otherf[comp]
  refsuba <- refsub
  refsuba[noadd] <- NA

  extrarefs <- refno[extra]

  # now, what if we should centre on the means?
  # there are two variants, either centre on the means in
  # both factors, or only in the one without a reference
  # we create a factor zfact which describes the groups
  # to mean.  

  
  # create a minimal environment for a function
  nm <- unlist(lapply(names(obj$fe),function(n) paste(n,levels(obj$fe[[n]]),sep='.')))
  fenv <- list(extrarefs=extrarefs,refsubs=refsubs,refsuba=refsuba,fef=fef,nm=nm)


  ef <- switch(as.character(opt),
         ln={
           local(function(v,addnames) {
             if(addnames) {
               names(v) <- nm
               attr(v,'extra') <- list(obs=obs,comp=comp,fe=fef,idx=idx)
             }
             v
           },list(obs=obs,comp=comp,fe=fef,idx=idx,nm=nm))
         },
         # one reference in each component
         ref={
           fenv$comp <- comp
           local(function(v,addnames) {
             esum <- sum(v[extrarefs])
             df <- v[refsubs]
             sub <- ifelse(is.na(df),0,df)
             df <- v[refsuba]
             add <- ifelse(is.na(df),0,df+esum)
             v <- v - sub + add
             if(addnames) {
               names(v) <- nm
               attr(v,'extra') <- list(obs=obs,comp=comp,fe=fef,idx=idx)
             }
             v
           },fenv)
         },
         zm={
           # now, what if we want zero-means on the other-factor?
           # we will then get an intercept for each component, and
           # zero means.  It's those which are in nosub, but partitioned
           # into components. We may do this faster now that we've
           # separated it from the ordinary 'ref'
           zfact <- comp
           zfact[!nosub] <- NA

           enames <- paste('icpt',1:ncomp,sep='.')
           zcomp <- factor(c(comp,1:ncomp))
           oo <- order(zcomp)
           fenv$oo <- oo
           fenv$zfact <- zfact
           fenv$zcomp <- zcomp[oo]
           fenv$enames <- enames
           fenv$obs <- c(obs,table(obj$cfactor))[oo]
           ef <- local(function(v,addnames) {
             esum <- sum(v[extrarefs])
             df <- v[refsubs]
             sub <- ifelse(is.na(df),0,df)
             df <- v[refsuba]
             add <- ifelse(is.na(df),0,df+esum)
             v <- v - sub + add
             means <- tapply(v,zfact,mean)
             mn <- means[zfact]
             mn <- ifelse(is.na(mn),0,mn)
             v <- v - mn
             v <- c(v,means)[oo]
             if(addnames) {
               names(v) <- c(nm,enames)[oo]
               attr(v,'extra') <- list(obs=obs,comp=zcomp)
             }
             v
           },fenv)
         },
         # one reference in each component, but zero-means in the other factor
         # and an intercept
         zm2={
           # both factors (but not the extra factors):
           # the interaction between comp and fef forms these groups,
           zfact <- interaction(comp,fef)
           # but skip the extra factors
           zfact[as.integer(comp) > ncomp] <- NA
           zfact <- factor(zfact)
           # and the means should be added to the intercepts
           # i.e. from the same component, add two and two
           # ifact should consist of the components of the
           # levels of zfact.
           ifact <- factor(as.integer(gsub('^([0-9]+).*','\\1',levels(zfact))),exclude=NA)
           enames <- paste('icpt',1:ncomp,sep='.')
           zcomp <- factor(c(comp,1:ncomp))
           oo <- order(zcomp)
           fenv$oo <- oo
           fenv$zfact <- zfact
           fenv$zcomp <- zcomp[oo]
           fenv$enames <- enames
           fenv$obs <- c(obs,table(obj$cfactor))[oo]
           fenv$ifact <- ifact
           ef <- local(function(v,addnames) {
             esum <- sum(v[extrarefs])
             df <- v[refsubs]
             sub <- ifelse(is.na(df),0,df)
             df <- v[refsuba]
             add <- ifelse(is.na(df),0,df+esum)
             v <- v - sub + add
             means <- tapply(v,zfact,mean)
             mn <- means[zfact]
             mn <- ifelse(is.na(mn),0,mn)
             v <- v - mn
             icpt <- tapply(means,ifact,sum)
             v <- c(v,icpt)[oo]
             if(addnames) {
               names(v) <- c(nm,enames)[oo]
               attr(v,'extra') <- list(obs=obs,comp=zcomp)
             }
             v
           },fenv)
         },

         stop(paste('estimable function',opt,'not recognized'))
         )

# try to byte compile the stuff
  o <- options(warn=-1)
  require('compiler',quietly=TRUE)
  options(o)
  if(exists('cmpfun'))
    ef <- cmpfun(ef,list(optimize=3))
  if(length(obj$fe) <= 2 && as.character(opt) != 'ln') 
    attr(ef,'verified') <- TRUE
  ef
} 


btrap <- function(alpha,obj,N=100,ef=NULL,eps=getOption('lfe.eps'),threads=getOption('lfe.threads')) {
  # bootstrap the stuff
  # bootstrapping is really to draw residuals over again, i.e. to change
  # the outcome.  Predictions of the estimated system are adjusted by
  # drawing from the residuals. We need a new r.h.s to replace
  # the current (I-P)(Y-Xbeta), i.e. (I-P)(Y-Xbeta+delta) where delta is resampled from
  # the residuals PY-PXbeta.  Then we adjust for degrees of freedom, which is component specific.


  if(is.null(ef))  {
    # use the one used with alpha
    ef <- attr(alpha,'ef')
  } else {
    if(is.character(ef)) ef <- efactory(obj,opt=ef)
    # redo the point estimates
    v <- ef(alpha[,'effect'],TRUE)
    alpha <- data.frame(effect=v)
    rownames(alpha) <- names(v)
    if(!is.null(attr(v,'extra'))) alpha <- cbind(alpha,attr(v,'extra'))
  }
  R <- obj$r.residuals-obj$residuals
  j <- 0 #avoid warning from check

  # We use the variance in PY-PXbeta to generate variation V (i.e. the full.residuals),
  # this is the outcome residuals, so the correct way to do it.
  # Then we compute W=(I-P)V.  Now, W is a vector which has constant
  # values for each level of each factor.  Couldn't we generate it by using a level
  # variation for the factors?  (This doesn't sound right, that's the coefficients...)
  # Then we adjust for the degrees of freedom, that's component specific
  # this may as well be done after kaczmarz solution, it makes no difference
  
  sefact <- sqrt((length(R)-1)/obj$df)
  sm <- mean(obj$residuals) # ought to be zero
  sefact <- 1  # do it at end instead
  if(length(obj$fe) > 1) {
    # need a factor into R assigning component number
    ac <- factor(unlist(lapply(obj$fe[1:2], function(f) obj$cfactor[match(levels(f),f)])))
    numpar <- table(ac)[ac]-1
    numobs <- table(obj$cfactor)[ac]
    o <- options(warn=-1)
    sefact <- as.vector(sqrt(length(R)/(numobs-numpar)))
    options(o)
  }

  # add extra factors
  if(length(obj$fe) > 2) {
    numex <- sum(unlist(lapply(obj$fe[3:length(obj$fe)],nlevels)))
    sefact <- c(sefact,rep(sqrt(length(R)/obj$df),numex))
  }
  # now, set those which we can't get standard error for to NA
  nose <- is.na(sefact)|is.infinite(sefact)
  sefact <- ifelse(nose,NA,sefact)
  smpdraw <- as.vector(obj$residuals)

  # Now, we want to do everything in parallel, so we should allocate up a set
  # of vectors, but we don't want to blow the memory.  Stick to allocating two
  # vectors per thread.  The threaded stuff can't be interrupted, so this is
  # an opportunity to control-c too.
  # hmm, up to 500 MB of vectors, we say, but no less than two per thread
  # (one per thread is bad for balance, if time to completion varies)
  maxB <- 500e6
  vpt <- max(2,as.integer(min(maxB/(length(R)*8),N)/threads))
  vpb <- vpt*threads
  blks <- as.integer(ceiling(N / vpb))
  newN <- blks*vpb
  vsum <- 0
  vsq <- 0
  start <- last <- as.integer(Sys.time())
  Rdup <- matrix(rep(R,vpb),length(smpdraw),vpb)
  gc()



  for(i in 1:blks) {
    rsamp <- sample(smpdraw,vpb*length(smpdraw),replace=TRUE)
    dim(rsamp) <- c(length(smpdraw),vpb)
    newR <- rsamp - demeanlist(rsamp,obj$fe,eps=eps,threads=threads) + Rdup
    v <- kaczmarz(obj$fe,newR,eps,threads=threads)*sefact
    efv <- apply(v,2,ef,addnames=FALSE)
    rm(v,rsamp,newR); gc()
    vsum <- vsum + rowSums(efv)
    vsq <- vsq + rowSums(efv**2)
    now <- as.integer(Sys.time())
    if(now-last > 300) {
      cat('...finished',i*vpb,'of',newN,'vectors in',now-start,'seconds\n')
      last <- now
    }
  }

  alpha[,'se'] <- sqrt(vsq/newN - (vsum/newN)**2)/(1-0.75/newN-7/32/newN**2-9/128/newN**3)

  return(alpha)
}


# check whether a function is estimable
is.estimable <- function(ef,fe,R=NULL,nowarn=FALSE) {
  if(!is.function(ef)) stop('ef must be a function')
  N <- sum(unlist(lapply(fe,nlevels)))
  if(is.null(R)) {
    # make a suitable residual
    nr <- length(fe[[1]])
    nc <- length(fe)
    vec <- unlist(lapply(fe,function(f) rnorm(nlevels(f))[f]))
    dim(vec) <- c(nr,nc)
    R <- rowSums(vec)
  }
  v1 <- ef(kaczmarz(fe,R,init=runif(N)),TRUE)
  v2 <- ef(kaczmarz(fe,R,init=runif(N)),TRUE)
  df <- max(abs(v1-v2))
  if(df > 1e-6) {
    bad <- which.max(abs(v1-v2))
    badname <- names(bad)
    if(!nowarn)
      warning('non-estimable function, largest error ',
              format(df,digits=1),' in coordinate ',bad, ' ("',badname,'")')
    return(FALSE)
  }
  TRUE
}
