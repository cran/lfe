# makematrix is a bit complicated. The purpose is to make model matrices for the various
# parts of the formulas.  The complications are due to the iv stuff.
# If there's an IV-part, its right hand side should be with the
# x. Their names are put in 'instruments'. Its left hand side goes in a separate entry 'ivy'



makematrix <- function(mf, contrasts=NULL, pf=parent.frame(),
                       clustervar=NULL, wildcard='n') {
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)

  # we should handle multiple lhs
  # but how?  model.frame() doesn't handle it, but we need
  # model.frame for subsetting and na.action, with the left hand side
  # included.  We create an artifical single lhs by summing the left hand
  # sides, just to get hold of the rhs.  Then we extract the left hand side

  # We need to remove the iv-spec from the Formula. It requires its own specification
  Form <- eval(mf[['formula']], pf)
  formenv <- environment(Form)
  Form <- as.Formula(Form)

  # If Form rhs is shorter than 4, extend it with zeros.
  # Then we avoid some special cases later
  numrhs <- length(Form)[2]
  
  # we can't just dot-update the iv-part, update will only keep the instruments
  if(numrhs < 2) Form <- update(Form, . ~ . | 0 | 0 | 0 | 0, drop=FALSE)
  else if(numrhs < 3) Form <- update(Form, . ~ . | . | 0 | 0 | 0 , drop=FALSE)
  else if(numrhs < 4) {
    # build from parts
    Form <- as.Formula(do.call(substitute, list(L ~ R1 | R2 | R3 | 0 | 0,
                                                list(L=formula(Form,lhs=NULL,rhs=0)[[2]],
                                                     R1=formula(Form,lhs=0,rhs=1)[[2]],
                                                     R2=formula(Form,lhs=0,rhs=2)[[2]],
                                                     R3=formula(Form,lhs=0,rhs=3)[[2]]))))
  } else if(numrhs < 5) {
    Form <- as.Formula(do.call(substitute, list(L ~ R1 | R2 | R3 | R4 | 0,
                                                list(L=formula(Form,lhs=NULL,rhs=0)[[2]],
                                                     R1=formula(Form,lhs=0,rhs=1)[[2]],
                                                     R2=formula(Form,lhs=0,rhs=2)[[2]],
                                                     R3=formula(Form,lhs=0,rhs=3)[[2]],
                                                     R4=formula(Form,lhs=0,rhs=4)[[2]]))))

  }


  if(numrhs > 5) stop("Formula can't have more than 5 parts")
# Make a suitable formula for a model frame. No tricky IV-spec
#  fullF <- formula(Form,lhs=NULL,rhs=0, drop=FALSE,collapse=TRUE,update=TRUE)
  fullF <- formula(Form,lhs=NULL,rhs=0, drop=FALSE)
  for(i in seq_len(length(Form)[2])) {
    f <- formula(Form,lhs=0,rhs=i,drop=FALSE)[[2]]
    if(i == 3) {
      if(identical(f,0)) next
      f <- as.Formula(f[[2]]) # skip '('
      f <- formula(f,collapse=TRUE, drop=FALSE)
      fullF <- update(fullF, formula(substitute(. ~ . + F1+F2, list(F1=f[[2]], F2=f[[3]]))), drop=FALSE, collapse=TRUE)
    } else {
      fullF <- update(fullF, formula(substitute(. ~ . + F, list(F=f))),drop=FALSE)
    }
  }

  usewild <- !identical(wildcard,'n') 
  dataenv <- new.env(parent=pf)
  if(usewild) {
    # we must evalaute the data argument, but we want to
    # avoid it being reevaluated when we eval(mf),
    # so put it in an environment.  We do it like this
    # to have a short name in mf[['data']] in case of errors.
    data <- eval(mf[['data']],pf)
    assign('..(@DATA@)..',data,dataenv)
    mf[['data']] <- as.name('..(@DATA@)..')
    wildnames <- colnames(data)
    rm(data)
    if(wildcard == 'R' || wildcard == 'G')
        wildnames <- unique(c(wildnames, rls(formenv)))
    rewild <- wildcard %in% c('r','R')
    fullF <- wildcard(fullF, wildnames, re=rewild)
  }
  environment(fullF) <- formenv

  mf[['formula']] <- fullF
  mf <- eval(mf, dataenv)
  rm(dataenv)
  naact <- na.action(mf)

#  if(is.null(mf$data)) data <- environment(mf[['formula']])
  # the factor list (rhs=2) needs special attention
  # it should be made into a model matrix, but treated specially.
  # It's a sum of terms like f + x:g

  fpart <- formula(Form, lhs=0, rhs=2)
  if(usewild) fpart <- wildcard(fpart,wildnames,re=rewild)
  ftm <- terms(fpart)

  # we make it into a call like
  # list(f=f, `x:g` = structure(g,x=x))
  # which we evaluate in the frame
  # make a function for ':'
  env <- new.env(parent = formenv)
# make  ':' a function of two arguments to do the interaction.
  assign(':', function(a,b) {
    anam <- deparse(substitute(a))
    bnam <- deparse(substitute(b))
    if(is.factor(a) && is.factor(b)) ret <- structure(interaction(a,b,drop=TRUE),xnam=bnam,fnam=anam)
    else if(is.factor(b)) ret <- structure(factor(b),x=a,xnam=anam,fnam=bnam)
    else if(is.factor(a)) ret <- structure(factor(a),x=b,xnam=bnam,fnam=anam)
    else stop('Error in term ',anam,':',bnam,'. Neither ',anam, ' nor ',bnam,' is a factor')
    ret
  }, env)

  fl <- lapply(attr(ftm,'term.labels'), function(tm) {
    f <- eval(parse(text=tm), mf, env)
    if(is.null(attr(f, 'fnam'))) factor(f) else f
  })
  
  names(fl) <- attr(ftm, 'term.labels')
  # Name the interactions with the matrix first, then the factor name
  names(fl) <- sapply(names(fl), function(n) {
    f <- fl[[n]]
    x <- attr(f,'x')
    if(is.null(x)) return(n)
    return(paste(attr(f,'xnam'),attr(f,'fnam'), sep=':'))
  })

  environment(Form) <- formenv
  if(is.null(clustervar)) {
    cluform <- terms(formula(Form, lhs=0, rhs=4))
    cluster <- lapply(eval(attr(cluform,'variables'), mf, pf), factor)
    names(cluster) <- attr(cluform,'term.labels')
    if(length(cluster) == 0) cluster <- NULL
  } else {
    # backwards compatible
    if(is.character(clustervar)) clustervar <- as.list(clustervar)
    if(!is.list(clustervar)) clustervar <- list(clustervar)
    cluster <- lapply(clustervar, function(cv) {
      if(!is.character(cv)) factor(cv) else factor(eval(as.name(cv),mf,formenv))
    })
  }

  ivform <- formula(Form,lhs=0, rhs=3, drop=FALSE)

  # Pick up IV instruments
  if(ivform[[1]] == as.name('~')) ivform <- ivform[[2]]
  if(ivform[[1]] == as.name('(')) ivform <- ivform[[2]]
  if(!identical(ivform,0)) {
    ivform <- as.Formula(ivform)
    if(length(ivform)[2] > 1) stop("Right hand side of IV-spec can't have multiple parts")
    inames <- as.character(attr(terms(formula(ivform, lhs=0, rhs=1)), 'variables'))[-1]
    environment(ivform) <- formenv
  } else {
    ivform <- NULL
    inames <- NULL
  }

  # then the fifth part, the controls
  form <- formula(Form, lhs=0, rhs=5, drop=TRUE)
  if(!identical(form[[2]],0)) {
    # always parse with intercept, remove it from matrix, so we never project out the intercept
    form <- formula(update(form, ~ . +1))
    if(usewild) form <- wildcard(form, wildnames, re=rewild)
    ctrlterms <- terms(form, data=mf)
    ctrl <- delete.icpt(model.matrix(ctrlterms, data=mf, contrasts.arg=contrasts))
    if(typeof(ctrl) != 'double') storage.mode(ctrl) <- 'double'
    if(ncol(ctrl) == 0) {
      ctrlnames <- ctrl <- NULL
    } else {
      ctrlnames <- colnames(ctrl)
    }
  } else {
    ctrl <- NULL
    ctrlnames <- NULL
  }

  # We have taken Form apart. Keep only exogenous variables
  Form <- formula(Form, lhs=NULL, rhs=1, drop=FALSE)
  environment(Form) <- formenv

# model.response doesn't work with multiple responses
#  y <- model.response(mf,"numeric") 

  form <- formula(Form, lhs=NULL, rhs=0, drop=FALSE)
  if(usewild) form <- wildcard(form, wildnames, re=rewild)
  y <- as.matrix(model.part(form, mf, lhs=NULL, rhs=0), rownames.force=FALSE)
  if(typeof(y) != 'double') storage.mode(y) <- 'double' 
  form <- formula(Form, lhs = 0, rhs = 1, collapse = c(FALSE, TRUE))
  if(usewild) form <- wildcard(form, wildnames, re=rewild)
  xterms <- terms(form, data=mf)
  x <- model.matrix(xterms, data=mf, contrasts.arg=contrasts)
  if(length(fl) > 0) {
    x <- delete.icpt(x)    
    icpt <- FALSE
  } else {
    icpt <- attr(xterms,'intercept') != 0
  }
  if(typeof(x) != 'double') storage.mode(x) <- 'double' 
  setdimnames(x, list(NULL, colnames(x)))

  if(!is.null(ivform)) {
    form <- formula(ivform, lhs=NULL, rhs=0, drop=FALSE)
    if(usewild) form <- wildcard(form, wildnames, re=rewild)
    ivy <-   as.matrix(model.part(form, mf, lhs=NULL, rhs=0), rownames.force=FALSE)
    if(typeof(ivy) != 'double') storage.mode(ivy) <- 'double' 
    form <- formula(ivform, lhs = 0, rhs = 1, collapse = c(FALSE, TRUE))
    if(usewild) form <- wildcard(form,wildnames,re=rewild)
    ivxterms <- terms(form, data=mf)
    # ivx should never contain an intercept
    ivx <- delete.icpt(model.matrix(ivxterms, data=mf, contrasts.arg=contrasts))
    if(typeof(ivx) != 'double') storage.mode(ivx) <- 'double' 
    setdimnames(ivx, list(NULL, colnames(ivx)))
  } else {
    ivy <- NULL
    ivx <- NULL
  }

  rm(mf) # save some memory

  # orig is necessary to compute the r.residuals, i.e. residuals without dummies
  # it's used in getfe() and btrap, but is of no use if we have ctrl variables
  TSS <- apply(y,2,var)*(nrow(y)-1)
  names(TSS) <- colnames(y)

  if(length(fl) != 0) {
    result <- edemeanlist(y=y, x=x, ivy=ivy, ivx=ivx, ctrl=ctrl, fl=fl)
    if(is.null(ctrl)) result$orig <- list(y=y, x=x, ivy=ivy, ivx=ivx)
  } else {
    result <- list(y=y, x=x, ivy=ivy, ivx=ivx, ctrl=ctrl)
  }
  rm(x,y,ivx,ivy,ctrl)  

  if(!is.null(result$ctrl)) {
    # pure control variables to project out
    # do ols, use the residuals as new variables
    y <- cbind(result$y,result$x,result$ivy,result$ivx)
    x <- result$ctrl
    result$ctrl <- NULL
#    fit <- .lm.fit(x,y)
# my own is much faster for large datasets
    fit <- newols(list(y=y,x=x), nostats=TRUE)
    resid <- as.matrix(fit$residuals)
    setdimnames(resid, list(NULL, colnames(y)))
    numctrl <- fit$rank
    rm(fit,x,y)

    result$y <- resid[,colnames(result$y), drop=FALSE]
    if(!is.null(result$x))   result$x <- resid[,colnames(result$x), drop=FALSE]
    if(!is.null(result$ivy)) result$ivy <- resid[,colnames(result$ivy), drop=FALSE]
    if(!is.null(result$ivx)) result$ivx <- resid[,colnames(result$ivx), drop=FALSE]
    rm(resid)
  } else {
    numctrl <- 0L
  }
  result$TSS <- TSS
  result$hasicpt <- icpt
  result$numctrl <- numctrl
  result$ctrlnames <- ctrlnames
  result$fl <- fl
  result$terms <- xterms
  result$cluster <- cluster
  result$formula <- Form
  result$ivform <- ivform
  result$inames <- inames
  result$na.action <- naact

  result
}


newols <- function(mm, stage1=NULL, pf=parent.frame(), nostats=FALSE, exactDOF=FALSE) {
  # We should project out the the factors from the y and x matrix, and do an ordinary lm.fit
#  if(length(mm$fl) != 0)
  if(!is.null(mm$orig))
      orig <- mm$orig
  else
      orig <- mm

  numctrl <- if(is.null(mm$numctrl)) 0 else mm$numctrl
  hasicpt <- if(is.null(mm$hasicpt)) FALSE else mm$hasicpt
  if(is.numeric(exactDOF)) {
    df <- exactDOF
    totvar <- nrow(mm$y) - df
  } else {
    # numrefs is also used later
    numrefs <- nrefs(mm$fl, compfactor(mm$fl), exactDOF) 
    totvar <- totalpvar(mm$fl)-numrefs + numctrl
    df <- nrow(mm$y)-totvar
  }

  # special case for no covariates

  if(is.null(mm$x) || ncol(mm$x) == 0) {
    
    z <- list(N=nrow(mm$x), r.residuals=orig$y,fe=mm$fl,p=totvar,Pp=0,cfactor=compfactor(mm$fl),
              na.action=mm$na.action, contrasts=mm$contrasts,
              fitted.values=orig$y - mm$y,
              df=df,
              nostats=FALSE,
              numctrl=numctrl,
              hasicpt=hasicpt,
              residuals=mm$y,clustervar=mm$cluster, call=match.call())
    z$df.residual <- z$df
    class(z) <- 'felm'
    return(z)
  }

  # lm.fit is an alternative.  Let's have a look at it later (didn't work before, but things have changed)
  # roll our own
  cp <- crossprod(mm$x)
  b <- crossprod(mm$x,mm$y)
  ch <- cholx(cp)
  badvars <- attr(ch,'badvars')

  z <- list()
  class(z) <- 'felm'
  if(is.null(badvars)) {
    beta <- backsolve(ch,backsolve(ch,b,transpose=TRUE))
    if(!nostats) z$inv <- chol2inv(ch)
  } else {
    beta <- matrix(NaN, nrow(cp), ncol(b))
    beta[-badvars,] <- backsolve(ch,backsolve(ch,b[-badvars,], transpose=TRUE))
    if(!nostats) {
      z$inv <- matrix(NA,nrow(cp),ncol(cp))
      z$inv[-badvars,-badvars] <- chol2inv(ch)
    }
  }
  rm(ch, b, cp)
  rownames(beta) <- colnames(orig$x)

  z$lhs <- colnames(beta) <- colnames(orig$y)
  z$hasicpt <- hasicpt
  z$TSS <- mm$TSS

  z$P.TSS <- apply(mm$y,2,var)*(nrow(mm$y)-1)
  names(z$P.TSS) <- colnames(mm$y)

  z$numctrl <- numctrl
  if(!nostats) inv <- nazero(z$inv)
  z$coefficients <- z$beta <- beta

  # what else is there to put into a felm object?
  z$Pp <- ncol(orig$x)
  z$N <- nrow(orig$x)
  z$p <- z$Pp - length(badvars) + numctrl

  nabeta <- nazero(beta)

  zfit <- mm$x %*% nabeta
  zresid <- mm$y - zfit
  
  z$response <- orig$y

  z$fitted.values <- mm$y - zresid
  z$residuals <- zresid
  z$contrasts <- mm$contrasts
  if(length(mm$fl) != 0) {
    z$r.residuals <- orig$y - orig$x %*% nabeta
  } else {
    z$r.residuals <- z$residuals
  }

  # the residuals should be the residuals from the original endogenous variables, not the predicted ones
  # the difference are the ivresid, which we must multiply by beta and subtract.
  # the residuals from the 2nd stage are in iv.residuals
  # hmm, what about the r.residuals?  We modify them as well. They are used in kaczmarz().

  if(!is.null(stage1)) {
    # we might need the centred response in condfstat()
    z$c.response <- mm$y
    fitnam <- makefitnames(stage1$lhs)
    ivresid <- stage1$residuals %*% nabeta[fitnam,,drop=FALSE]
    z$iv.residuals <- z$residuals
    z$residuals <- z$residuals - ivresid    
    z$r.iv.residuals <- z$r.residuals
    z$r.residuals <- z$r.residuals - ivresid
    z$endovars <- stage1$lhs
  }

  z$terms <- mm$terms
  z$cfactor <- compfactor(mm$fl)
  totlev <- totalpvar(mm$fl)

  if(is.numeric(exactDOF)) {
    z$df <- exactDOF
    numdum <- z$N - z$p - z$df
    z$numrefs <- totlev - numdum
  } else {
    numdum <- totlev - numrefs
    z$numrefs <- numrefs
    z$df <- z$N - z$p - numdum
  }
  z$df.residual <- z$df
  z$rank <- z$N - z$df
  z$exactDOF <- exactDOF

  z$fe <- mm$fl
# should we subtract 1 for an intercept?
# a similar adjustment is done in summary.felm when computing rdf
  z$p <- z$p + numdum #- 1  
  z$xp <- z$p
  z$na.action <- mm$na.action
  class(z) <- 'felm'
  cluster <- mm$cluster
  z$clustervar <- cluster
  z$stage1 <- stage1
  if(nostats) {
    z$nostats <- TRUE
    return(z)
  }

  z$nostats <- FALSE

# then we go about creating the covariance matrices and tests
# if there is a single lhs, they are just stored as matrices etc
# in z.  If there are multiple lhs, these quantities are inserted
# in a list z$STATS indexed by z$lhs
# indexed by the name of the lhs

  vcvnames <- list(rownames(beta), rownames(beta))
  Ncoef <- nrow(beta)

  singlelhs <- length(z$lhs) == 1
  # preallocate STATS 
  if(!singlelhs) z$STATS <- list()
  z$STATS <- list()

  for(lhs in z$lhs) {
    res <- z$residuals[,lhs]

# when multiple lhs, vcvfactor is a vector
# we need a list of vcvs in this case

    vcv <- sum(res**2)/z$df * z$inv
    setdimnames(vcv, vcvnames)
      z$STATS[[lhs]]$vcv <- vcv
    if(singlelhs) z$vcv <- vcv

  # We should make the robust covariance matrix too.
  # it's inv * sum (X_i' u_i u_i' X_i) * inv
  # where u_i are the (full) residuals (Wooldridge, 10.5.4 (10.59))
  # i.e. inv * sum(u_i^2 X_i' X_i) * inv
  # for large datasets the sum is probably best computed by a series of scaled
  # rank k updates, i.e. the dsyrk blas routine, we make an R-version of it.
  # need to check this computation, the SE's are slightly numerically different from Stata's.
  # it seems stata does not do the small-sample adjustment
    dfadj <- z$N/z$df

  # Now, here's an optimzation for very large xz. If we use the wcrossprod and ccrossprod
  # functions, we can't get rid of xz, we end up with a copy of it which blows away memory.
  # we need to scale xz with the residuals in xz, but we don't want to expand res to a full matrix,
  # and even get a copy in the result.
  # Thus we modify it in place with a .Call. The scaled variant is also used in the cluster computation.

    xz <- mm$x
    rscale <- ifelse(res==0,1e-40,res)  # make sure nothing is zero
    # This one scales the columns without copying
    # For xz, remember to scale it back, because we scale directly into
    # mm$x
    .Call(C_scalecols, xz, rscale)
    # compute inv %*% crossprod(xz) %*% inv
    # via a blas dsyrk. Save some memory 
    meat <- matrix(0, Ncoef, Ncoef)
    .Call(C_dsyrk,0.0,meat,dfadj,xz)
    rvcv <- inv %*% meat %*% inv
    setdimnames(rvcv, vcvnames)

    z$STATS[[lhs]]$robustvcv <- rvcv
    if(singlelhs) z$robustvcv <- rvcv

    rm(meat, rvcv)

  # then the clustered covariance matrix
    if(!is.null(cluster)) {
      method <- attr(cluster,'method')
      if(is.null(method)) method <- 'cgm'
      dfadj <- (z$N-1)/z$df
      d <- length(cluster)
      if(method == 'cgm' || TRUE) {
        meat <- matrix(0,Ncoef,Ncoef)
        for(i in 1:(2^d-1)) {
          # Find out which ones to interact
          iac <- as.logical(intToBits(i))[1:d]
          # odd number is positive, even is negative
          sgn <- 2*(sum(iac) %% 2) - 1
          # interact the factors
          ia <- factor(do.call(paste,c(cluster[iac],sep='\004')))
          adj <- sgn*dfadj*nlevels(ia)/(nlevels(ia)-1)
          .Call(C_dsyrk,1.0,meat,adj,rowsum(xz,ia))
        }
        cvcv <- inv %*% meat %*% inv
        setdimnames(cvcv, vcvnames)
        z$STATS[[lhs]]$clustervcv <- cvcv
        if(singlelhs) z$clustervcv <- cvcv
        rm(meat,cvcv)

        ## } else if(method == 'gaure') {
        ##   .Call(C_scalecols, xz, 1/rscale)
        ##   meat <- matrix(0,nrow(z$vcv),ncol(z$vcv))
        ##   dm.res <- demeanlist(res,cluster)
        ##   skel <- lapply(cluster, function(f) rep(0,nlevels(f)))
        ##   means <- relist(kaczmarz(cluster,res-dm.res), skel)
        ##   scale <- ifelse(dm.res==0,1e-40, dm.res)
        ##   .Call(C_scalecols, xz, scale)
        ##   .Call(C_dsyrk, 1, meat, dfadj, xz)
        ##   .Call(C_scalecols, xz, 1/scale)
        ##   for(i in seq_along(cluster)) {
        ##     rs <- rowsum(xz, cluster[[i]])
        ##     adj <- nlevels(cluster[[i]])/(nlevels(cluster[[i]])-1)
        ##     .Call(C_scalecols, rs, means[[i]])
        ##     .Call(C_dsyrk, 1, meat, dfadj*adj, rs)
        ##   }
        ##   rm(xz,rs)
        ##   z$clustervcv <- inv %*% meat %*% inv
        ##   rm(meat)
      } else {
        stop('unknown multi way cluster algorithm:',method)
      }
      
      
      z$STATS[[lhs]]$cse <- sqrt(diag(z$STATS[[lhs]]$clustervcv))
      z$STATS[[lhs]]$ctval <- z$coefficients[,lhs]/z$STATS[[lhs]]$cse
      z$STATS[[lhs]]$cpval <- 2*pt(abs(z$STATS[[lhs]]$ctval),z$df,lower.tail=FALSE)
      
      if(singlelhs) {
        z$cse <- z$STATS[[lhs]]$cse
        z$ctval <- z$STATS[[lhs]]$ctval
        z$cpval <- z$STATS[[lhs]]$cpval
      } 
    }
    z$STATS[[lhs]]$se <- sqrt(diag(z$STATS[[lhs]]$vcv))
    z$STATS[[lhs]]$tval <- z$coefficients[,lhs]/z$STATS[[lhs]]$se
    z$STATS[[lhs]]$pval <- 2*pt(abs(z$STATS[[lhs]]$tval),z$df,lower.tail=FALSE)
    
    z$STATS[[lhs]]$rse <- sqrt(diag(z$STATS[[lhs]]$robustvcv))
    z$STATS[[lhs]]$rtval <- z$coefficients[,lhs]/z$STATS[[lhs]]$rse
    z$STATS[[lhs]]$rpval <- 2*pt(abs(z$STATS[[lhs]]$rtval),z$df,lower.tail=FALSE)
    
    if(singlelhs) {
      z$se <- z$STATS[[lhs]]$se
      z$tval <- z$STATS[[lhs]]$tval
      z$pval <- z$STATS[[lhs]]$pval

      z$rse <- z$STATS[[lhs]]$rse
      z$rtval <- z$STATS[[lhs]]$rtval
      z$rpval <- z$STATS[[lhs]]$rpval
    }

    # reset this for next lhs
    .Call(C_scalecols, xz, 1/rscale)
  }
  z
}

felm <- function(formula, data, exactDOF=FALSE, subset, na.action, contrasts=NULL,...) {

  knownargs <- c('iv', 'clustervar', 'cmethod', 'keepX', 'nostats', 'wildcard')
  keepX <- FALSE
  cmethod <- 'cgm'
  iv <- NULL
  clustervar <- NULL
  nostats <- FALSE
  wildcard <- 'n'
  deprec <- c('iv', 'clustervar')
  mf <- match.call(expand.dots = FALSE)

  # Currently there shouldn't be any ... arguments
  # check that the list is empty

#  if(length(mf[['...']]) > 0) stop('unknown argument ',mf['...'])
  
  # When moved to the ... list, we use this:
  # we do it right away, iv and clustervar can't possibly end up in ... yet, not with normal users


  args <- list(...)
  ka <- knownargs[pmatch(names(args),knownargs, duplicates.ok=FALSE)]
  names(args)[!is.na(ka)] <- ka[!is.na(ka)]
  dpr <- deprec[match(ka, deprec)]
  if(any(!is.na(dpr))) {
    bad <- dpr[which(!is.na(dpr))]
    warning('Argument(s) ',paste(bad,collapse=','), ' are deprecated and will be removed, use multipart formula instead')
  }
  env <- environment()
  lapply(intersect(knownargs,ka), function(arg) assign(arg,args[[arg]], pos=env))

  if(!(cmethod %in% c('cgm','gaure'))) stop('Unknown cmethod: ',cmethod)

  # also implement a check for unknown arguments
  unk <- setdiff(names(args), knownargs)
  if(length(unk) > 0) stop('unknown arguments ',paste(unk, collapse=' '))

  # backwards compatible
  Gtm <- terms(formula(as.Formula(formula), rhs=1), specials='G')
  if(!is.null(attr(Gtm,'specials')$G) || !is.null(iv)) {
    mf <- match.call(expand.dots=TRUE)
    mf[[1L]] <- quote(..oldfelm)
    return(eval.parent(mf))
  }

  mm <- makematrix(mf, contrasts, pf=parent.frame(), clustervar, wildcard=wildcard)
  ivform <- mm$ivform

  if(is.null(ivform)) {
    # no iv, just do the thing
    z <- newols(mm, nostats=nostats[1], exactDOF=exactDOF)
    if(keepX) z$X <- if(is.null(mm$orig)) mm$x else mm$orig$x
    z$call <- match.call()
    return(z)
  }


  ########### Instrumental variables ############
  # pick up the full formula
  fullform <- mm$formula

  if(length(nostats) == 2)
      nost1 <- nostats[2]
  else
      nost1 <- nostats[1]


  # Now, we must build a model matrix with the endogenous variables
  # on the left hand side, the exogenous and instruments on the rhs.
  # we have already centred everything in mm. However, we must
  # rearrange it.
  # in the first stage we should have the iv left hand side on
  # the lhs, the exogenous and instruments on the rhs.
  # mm$x is an ok rhs. The y must be replaced by the ivy
  
  ivars <- colnames(mm$ivx)
  exlist <- colnames(mm$x)

  mm1 <- mm[names(mm) %in% c('fl','terms','cluster', 'numctrl',
                             'hasicpt','na.action','contrasts')]
  mm1$y <- mm$ivy
  mm1$x <- cbind(mm$x, mm$ivx)
  mm1$orig$y <- mm$orig$ivy
  mm1$orig$x <- cbind(mm$orig$x, mm$orig$ivx)

  fitnames <- makefitnames(colnames(mm$ivy))

  z1 <- newols(mm1, nostats=nost1, exactDOF=exactDOF)
  rm(mm1)

  if(!nost1) {
    z1$iv1fstat <- lapply(z1$lhs, function(lh) waldtest(z1,ivars, lhs=lh))
    names(z1$iv1fstat) <- z1$lhs
    z1$rob.iv1fstat <- lapply(z1$lhs, function(lh) waldtest(z1,ivars,type='robust', lhs=lh))
    names(z1$rob.iv1fstat) <- z1$lhs
  }

  z1$instruments <- ivars
  z1$centred.exo <- mm$x
  z1$ivx <- mm$ivx 
  z1$ivy <- mm$ivy

  # then second stage.  This is a bit more manipulation
  # We must make an mm2 with the original lhs, and with
  # the exogenous variables plus the the predicted endogenous from z1 on the rhs
  # we must set the names of the exogenous variables

  mm2 <- mm[names(mm) %in% c('fl','terms','cluster','numctrl','hasicpt',
                             'na.action','contrasts', 'TSS')]
  mm2$x <- cbind(mm$x, z1$fitted.values)
  setdimnames(mm2$x, list(NULL, c(exlist,fitnames)))
  mm2$y <- mm$y
  if(!is.null(mm$orig)) {
    mm2$orig <- list(x=cbind(mm$orig$x, z1$fitted.values), y=mm$orig$y)
    setdimnames(mm2$orig$x, list(NULL, c(exlist,fitnames)))
  }
  rm(mm)  # save some memory
  z2 <- newols(mm2, stage1=z1, nostats=nostats[1], exactDOF=exactDOF)
  if(keepX) z2$X <- if(is.null(mm2$orig)) mm2$x else mm2$orig$x
  rm(mm2)

  z2$call <- match.call()
  z2
}
