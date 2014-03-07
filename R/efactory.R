efactory <- function(obj, opt='ref', ...) {

  # only factors without covariates are relevant to analyze
  purefes <- sapply(obj$fe, function(f) is.null(attr(f,'x')))
  pfe <- obj$fe[purefes]

#  allnm <- unlist(lapply(names(obj$fe),function(n) paste(n,levels(obj$fe[[n]]),sep='\003')))
  allnm <- unlist(lapply(names(obj$fe),function(n) xlevels(n,obj$fe[[n]],sep='\003')))
  
# the names of the dummies, e.g. id.4 firm.23
#  nm <- unlist(lapply(names(pfe),function(n) paste(n,levels(pfe[[n]]),sep='\003')))
  nm <- unlist(lapply(names(pfe),function(n) xlevels(n,pfe[[n]],sep='\003')))
  # create an index where the pure fe's belong in the full array
  allpos <- match(nm,allnm)
  mkallvec <- function(x) {res <- rep(NA,length(allnm)); res[allpos] <- allpos[x]; res;}
# how many obervations for each level
  lobs <- lapply(pfe,table)
  obs <- unlist(lobs)  
  names(obs) <- unlist(lapply(names(lobs), function(n) paste(n,'\003',names(lobs[[n]]),sep='')))

#  allobs <- unlist(lapply(obj$fe,table))
  allobs <- unlist(lapply(obj$fe,function(f) {
    x <- attr(f,'x')
    if(is.null(x)) return(table(f))
    if(!is.matrix(x)) return(table(f))
    return(rep(table(f), ncol(x)))
  }))  

  if(length(pfe) == 2) {
    # now, find the component of each parameter, i.e. each level.  We do this
    # by finding the first occurence of each level, i.e. match(levels(f),f)
    comp <- factor(unlist(lapply(pfe, function(f) obj$cfactor[match(levels(f),f)])))
    ncomp <- nlevels(comp)
  } else if(length(pfe) > 2) {
    # we should formally assign unique component numbers for factors beyond the second
    comp <- factor(unlist(lapply(pfe[1:2], function(f) obj$cfactor[match(levels(f),f)])))
    ncomp <- nlevels(comp)
    exlvls <- (nlevels(comp)+1):(nlevels(comp)+1 + length(pfe)-3)
    comp <- as.factor(c(comp,unlist(mapply(rep,exlvls,unlist(lapply(pfe[3:length(pfe)],nlevels))))))
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

  fef <- factor(unlist(lapply(names(pfe),function(n) rep(n,nlevels(pfe[[n]])))))
#  allfef <- factor(unlist(lapply(names(obj$fe),function(n) rep(n,nlevels(obj$fe[[n]])))))
  allfef <- factor(unlist(lapply(names(obj$fe),function(n) nxlevels(n,obj$fe[[n]]))))

  # level of the factor
#  idx <- factor(unlist(lapply(obj$fe,function(f) levels(f))))
  idx <- factor(unlist(lapply(obj$fe,function(f) {
    x <- attr(f,'x')
    if(is.null(x) || !is.matrix(x)) return(levels(f))
    return(rep(levels(f), ncol(x)))
  })))
  # then figure out in which factor the reference is
  # make sure to allow '.' in factor names
  rf <- sub('(^.*)\003..*$','\\1',refnames)
  # now, create a refsubs which is the ones to be subtracted
  # each refsub belonging to somthing else than the reference factor
  # should be NA'ed.
  if(length(pfe) > 2) {
    extra <- (length(refno)-length(pfe)+3):length(refno)
    sw <- c(names(pfe)[c(2,1)],rep('.NA',length(pfe)-2))
  } else {
    swap <- if(length(pfe) == 2) c(2,1) else 1
    sw <- names(pfe)[swap]
    extra <- integer(0)
  }
  names(sw) <- names(pfe)
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

  extrarefs <- allpos[extrarefs]
  refsubs <- mkallvec(refsubs)
  refsuba <- mkallvec(refsuba)
  obs <- allobs
  fef <- allfef
#  nm <- unlist(lapply(names(obj$fe),function(n) paste(n,levels(obj$fe[[n]]),sep='.')))
  nm <- unlist(lapply(names(obj$fe),function(n) xlevels(n,obj$fe[[n]],sep='.')))
#  allcomp <- rep(0,sum(sapply(obj$fe,nlevels)))
  allcomp <- rep(0,length(allnm))
  allcomp[allpos] <- comp
  comp <- allcomp
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
           fenv$mkallvec <- mkallvec
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
  if(length(pfe) <= 2 && as.character(opt) != 'ln' && all(purefes)) 
    attr(ef,'verified') <- TRUE
  ef
} 
