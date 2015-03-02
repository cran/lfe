ivbootstrap <- function(z, x, y, quantiles=0.95, N=100L, cluster=NULL) {
  # estimate bias as E((z' x)^{-1} z' eps)
  N <- max(N)
  if(!is.null(cluster)) {
    if(length(cluster) > 1)
        warning('Only a single cluster is supported for IV bootstrap, using ',
                names(cluster)[[1]])
    clu <- cluster[[1]]
    iclu <- as.integer(clu)
  }
#  zortho <- orthonormalize(z)

  pint <- getOption('lfe.pint')
  start <- last <- Sys.time()
  # hmm, can we reduce the number of instruments?

  n <- 0
  bias <- replicate(N,{
    n <<- n + 1
    now <- Sys.time()
    if(is.numeric(pint) && pint > 0 && now-last > pint) {
      message(date(), ' Iteration ', n , ' of ', N , ' in IV bootstrap')
      last <<- now
    }
    if(is.null(cluster)) {
      # resampling observations for indep residuals
      s <- sample(nrow(z),replace=TRUE)
    } else {
      # resample entire levels
      cl <- sort(sample(nlevels(clu), replace=TRUE))
      # find a faster way to do this:
      # s <- sort(unlist(sapply(cl, function(ll) which(clu==ll))))
      s <- NULL
      while(length(cl) > 0) {
        s <- c(s,which(iclu %in% cl))
        cl <- cl[c(1L,diff(cl)) == 0]
      }
    }
    
    # draw new instruments
    zortho <- orthonormalize(z[s,,drop=FALSE])
    # draw new X, and project it
    zX <- crossprod(zortho, x[s,,drop=FALSE])
    tryCatch(solve(crossprod(zX),
                   crossprod(zX, crossprod(zortho, y[s,,drop=FALSE]))),
             error=function(e) {warning(e);NULL})
  })
  if(start != last) cat('\n')

  if(is.list(bias)) {
    # some of them returned NULL, so replicate returned a list
    # discard those NULLs
    bias <- simplify2array(bias[!sapply(bias,is.null)])
    N <- dim(bias)[3]
  }
  # estimate quantiles of the bias

  if(is.null(quantiles)) {
    res <- bias
    res <- aperm(res, c(1,(3:length(dim(res))),2))
  } else {
    qname <- paste(100*round(quantiles,3), '%',sep='')
    dm <- 1:(length(dim(bias))-1)
    res <- apply(bias,dm,function(s) quantile(s,probs=quantiles,na.rm=TRUE,type=4))
    if(length(quantiles) == 1) {
      dmn <- dimnames(res)
      dim(res) <- c(1,dim(res))
      dimnames(res) <- c(list(paste(round(quantiles,3),'%',sep='')),dmn)
    }
    res <- aperm(res,c(2,1,3:length(dim(res))))
  }

  structure(res, q=quantiles, samples=N)
}

  # From "A weak instrument F-test in linear iv models with multiple
  # " endogenous variables", Sanderson & Windmeijer, Disc. Paper 14/644 U of Bristol, 2014
  # pp 22-23.  There's an error in how tilde delta is computed at p. 23.
  # it should be \hat X_{-j}, not X_{-j}.  I.e. coefs for \hat X_{-j} when regressing
  # x_j.  But we should predict with with X_{-j}
  # I.e. estimate delta in x_j = \hat X_{-j} * delta + eps
  # Regress x_j - X_{-j}*delta = Z * kappa + xi
  # wald test on kappa, with df = max(1,kz-#endog+1)
condfstat <- function(object, type='default', quantiles=0.0, bN=100L) {
  est <- object
  st1 <- est$stage1
  if(is.null(st1))
      stop('Conditional F statistic only makes sense for iv-estimation')


  if(is.null(type)) {
    types <- c('iid','robust')
    if(!is.null(est$clustervar)) types <- c(types,'cluster')
  } else {
    if(identical(type,'default'))
        types <- if(is.null(est$clustervar)) 'iid' else 'cluster'
    else
        types <- type
  }

  if(length(st1$lhs) == 1) {
    # only a single endogenous variable
    # reduce to ordinary F-test
    df1 <- nrow(st1$coefficients)
    result <- as.matrix(sapply(types, function(typ) {
       waldtest(st1,st1$instruments, df1=df1, type=typ)['F']
    }))
    dimnames(result) <- list(st1$lhs,paste(types,'F'))
    return(structure(t(result),df1=df1))
  }


  # first, transform away all the exogenous variables from
  # instruments, endogenous variables and predicted endogenous variables
  # there may be an intercept in ivx, we remove it
  keep <- !(colnames(st1$ivx) %in% '(Intercept)')
  if(all(keep)) ivx <- st1$ivx else ivx <- st1$ivx[,keep, drop=FALSE]
  inames <- colnames(ivx)
  y <- cbind(st1$ivy, ivx, st1$fitted.values, est$c.response)
  fitnames <- makefitnames(colnames(st1$fitted.values))
  setdimnames(y, list(NULL, c(colnames(st1$ivy), inames, fitnames, colnames(est$c.response))))
  mm <- list(y=y, x=st1$centred.exo)
  tvars <- newols(mm, nostats=TRUE)$residuals
  rm(y,mm)
  # should we estimate the relative bias?
  if(any(quantiles > 0) || is.null(quantiles)) {
    # use the predicted variables as instruments?
    z <- tvars[,colnames(ivx),drop=FALSE]
#    z <- xhat
    x <- tvars[,colnames(st1$ivy),drop=FALSE]
    y <- tvars[,colnames(est$c.response), drop=FALSE]
    bias <- ivbootstrap(z,x,y,
                   quantiles=quantiles,N=bN,cluster=est$clustervar)
    rm(z,x,y)

    # Now, bias contains the bias distribution
    # we subtract it from the estimate in object$coefficients
    quant <- bias
#    cf <- object$coefficients[fitnames,,drop=FALSE]
#    quant <- sapply(seq_len(ncol(cf)), function(i) cf[,i] + bias[,,i,drop=FALSE])
#    attributes(quant) <- attributes(bias)
            
    quant <- drop(quant)
  } else {
    quant <- NULL
  }

  resid <- sapply(st1$lhs, function(ev) {
    # do an ols on each
    evrest <- st1$lhs[!(st1$lhs %in% ev)]
    hatev <- makefitnames(ev)
    hatevrest <- makefitnames(evrest)
    Xlhs <- tvars[,ev,drop=FALSE]
    Xhatrest <- tvars[,hatevrest,drop=FALSE]
    Xrest <- tvars[,evrest,drop=FALSE]
    Xlhs - Xrest %*% solve(crossprod(Xhatrest), t(Xhatrest) %*% Xlhs)
  })

  setdimnames(resid, list(NULL, st1$lhs))


  # then regress on the instruments
  mm <- list(y = resid, x=tvars[,inames], cluster=est$clustervar)
  z <- newols(mm, nostats=FALSE)

  df1 <- nrow(z$coefficients)-length(z$lhs)+1
  result <- as.matrix(sapply(types, function(typ) {
    sapply(z$lhs, function(lh) waldtest(z,inames,lhs=lh, df1=df1, type=typ)['F'])
  }))
  dimnames(result) <- list(z$lhs,paste(types,'F'))
  structure(t(result),df1=df1, quantiles=quant)
}
 
oldcondfstat <- function(object, type='default') {
  est <- object
  if(is.null(est$stage1))
      stop('Conditional F statistic only makes sense for iv-estimation')

  # for each endogenous variable x_j, we move it to the lhs, and estimate the
  # residuals with the other endogenous vars as explanatory variables.

  # We put each of these residuals on the lhs in the 1st stage, and estimate
  # the effects of the instruments.  We do a Wald test on these estimates to obtain
  # a conditional test for each endogenous variable.

  st1 <- est$stage1
  if(length(st1$lhs) == 1) {
    # only a single endogenous variable
    # reduce to ordinary F-test
    W <- waldtest(st1,st1$instruments)
    return(structure(W['F'],df1=W['df1']))
  }

  cl1 <- st1$call

  newcl <- cl <- est$call
  Form <- as.Formula(cl[['formula']])
  baseform <- as.Formula(formula(Form,lhs=0,rhs=1:2))
  ivForm <- as.Formula(formula(Form,lhs=0,rhs=3)[[2]][[2]])
  endogForm <- formula(ivForm, lhs=NULL, rhs=0)[[2]]
  endogvars <- all.vars(endogForm)
  instrvars <- st1$instruments    #all.vars(formula(ivForm, lhs=0, rhs=NULL)[[2]])

  fitvars <- st1$endogvars

  # Now, the residuals should be put into the lhs of the 1st stage
  # The 1st stage formula is already ok, but we should make sure
  # the lhs is picked up from the residuals we create.
  # We create a new environment based on the old, but our new
  # residuals into the new environment

  fitenv <- environment(est$st2call[['formula']])
  # an environment to store the new residual variables

  # we can get rid of the exogeneous residuals by
  # replacing the instruments and endogeneous variables
  # by their residuals when regressing with the exogeneous variables
  # If K is the number of endogenous variables, we do K+1 centerings
  # of the exogenous variables below.  If we project them out first,
  # we will just do 1 centering.
  # Collect the endogenous variables and instruments and their projections on the rhs
  # Only the exogenous variables on the rhs. Store residuals in noexo environment.

  noexo <- new.env()

  # put endogenous variables, their predictions from 1st stage, and the
  # instruments on the lhs, with the exogenous variables on the rhs
  cvars <- c(endogvars, fitvars, instrvars)
  cForm <- as.Formula(formula(paste(paste(cvars, collapse='|'),'~0')))

  cForm <- as.Formula(update(cForm, as.formula(substitute(. ~ X,
                                                          list(X=baseform[[2]])))))
  projcall <- cl
  environment(cForm) <- fitenv
  projcall[['formula']] <- cForm
  projcall[['nostats']] <- TRUE
  projest <- eval(projcall, envir=est$parent.frame)

  # store residuals in environment under name 'x..(noexo)..
  # We store the residuals we need for the next regression
  # i.e. the x_j - X_{-j} \delta, under the name
  # of the endogenous variable.  That's all we need.
  # and the projected instruments.

  lhsvars <- colnames(st1$residuals)
  endogfitvars <- paste(lhsvars,'(fit)',sep='')
  for(ev in seq_along(lhsvars)) {
    # do an OLS to compute delta, but implement just with crossprod on the right matrix
    evnam <- lhsvars[ev]
    restfitnam <- endogfitvars[-ev]
    restnam <- lhsvars[-ev]
    orig <-  projest$residuals[,evnam,drop=FALSE]
    Xjfit <- projest$residuals[,restfitnam,drop=FALSE]
    Xj <- projest$residuals[,restnam,drop=FALSE]
    delta <- solve(crossprod(Xjfit), t(Xjfit) %*% orig)
    resid <- orig - Xj %*% delta
    # these are the left hand sides in the stage 2 regression below
    resnam <- paste(evnam,'..(noexo)..',sep='')
    assign(resnam, resid, envir=noexo)
  }
  # store the residual instruments in noexo also

  for(instr in instrvars) {
    resnam <- paste(instr,'..(noexo)..',sep='')
    assign(resnam, projest$residuals[,instr,drop=FALSE], envir=noexo)
  }
  # ditch projest to save memory, we don't need it anymore
  rm(projest)

  # compute concentration parameter matrix mu
  # From Stock, Wright, Yogo section 4.2 p .522
  # do this when I have the time to do it.
  
  # Then do the regression on the endog. residuals in resenv
  # on the lhs, and projected instruments from
  # projest on the rhs
  # make a formula
  # Should we do clustering?
  resendog <- paste('`',lhsvars,'..(noexo)..`',sep='')
  resinst <- paste('`',instrvars,'..(noexo)..`',sep='')
  if(length(Form)[2] == 4) {
    cluform <- formula(formula(Form,lhs=0,rhs=4))[[2]]
    F1 <- as.Formula(formula(paste(paste(resendog,collapse='|'),'~',
                                   paste(resinst,collapse='+'),'+0|0|0|0')))
    newform <- update(F1,as.formula(substitute(. ~ . |.|.|C, list(C=cluform))))

  } else {
    newform <- as.Formula(formula(paste(paste(resendog,collapse='|'),'~',
                                      paste(resinst,collapse='+'),'+0')))
  }

  # everything except the cluster var is in noexo
  # the cluster var is still in the original data frame, or in the parent
  # of noexo. We must call felm with the original data. Hopefully our
  # constructed names ..(residual).. are not present there
  environment(newform) <- noexo
  m <- match(c('formula','data'), names(cl),0L)
  newcl <- cl[c(1L,m)]
  newcl[['formula']] <- newform
  e1 <- eval(newcl,est$parent.frame)

  df1 <- length(instrvars)-length(lhsvars)+1
  res <- sapply(e1$lhs, function(lh) {
    waldtest(e1,resinst,lhs=lh, df1=df1, type=type)['F']
  })
  attr(res,'df1') <- df1
  names(res) <- lhsvars
  res
}
