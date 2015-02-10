
btrap <- function(alpha,obj,N=100,ef=NULL,eps=getOption('lfe.eps'),
                  threads=getOption('lfe.threads'), robust=FALSE,
                  cluster=NULL, lhs=NULL) {
  # bootstrap the stuff. The name 'btrap' is chosen to instill a feeling of being trapped
  # (in some long running stuff which will never complete)
  # bootstrapping is really to draw residuals over again, i.e. to change
  # the outcome.  Predictions of the estimated system are adjusted by
  # drawing from the residuals. We need a new r.h.s to replace
  # the current (I-P)(Y-Xbeta), i.e. (I-P)(Y-Xbeta+delta) where delta is resampled from
  # the residuals PY-PXbeta.  Then we adjust for degrees of freedom, which is component specific.
    
  if(is.logical(cluster)) {
      if(cluster) cluster <- obj$clustervar else cluster <- NULL
  } else if(!is.null(cluster)) {
      if(is.list(cluster))
          cluster <- lapply(cluster,factor)
      else
          cluster <- list(factor(cluster))
  }

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

  if(is.null(lhs))
      R <- obj$r.residuals-obj$residuals
  else
      R <- obj$r.residuals[,lhs]-obj$residuals[,lhs]
  
  # We use the variance in PY-PXbeta to generate variation V (i.e. the full.residuals),
  # this is the outcome residuals, so the correct way to do it.
  # Then we compute W=(I-P)V.  Now, W is a vector which has constant
  # values for each level of each factor.  Couldn't we generate it by using a level
  # variation for the factors?  (This doesn't sound right, that's the coefficients...)
  # Then we adjust for the degrees of freedom, that's component specific
  # this may as well be done after kaczmarz solution, it makes no difference
  Rvec <- as.vector(R)  
  sefact <- sqrt(obj$N/obj$df)
  if(is.null(lhs)) {
    smpdraw <- as.vector(obj$residuals)
  }  else {
    smpdraw <- as.vector(obj$residuals[,lhs])
  }


  # Now, we want to do everything in parallel, so we should allocate up a set
  # of vectors, but we don't want to blow the memory.  Stick to allocating two
  # vectors per thread.  The threaded stuff can't be interrupted, so this is
  # an opportunity to control-c too.
  # hmm, up to 500 MB of vectors, we say, but no less than two per thread
  # (one per thread is bad for balance, if time to completion varies)
  # divide by two because we use a copy in the demeanlist step.
  maxB <- getOption('lfe.bootmem')*1e6/2
  vpt <- max(2,as.integer(min(maxB/(length(R)*8),N)/threads))
  vpb <- vpt*threads
  blks <- as.integer(ceiling(N / vpb))
  newN <- blks*vpb
  vsum <- 0
  vsq <- 0
  start <- last <- as.integer(Sys.time())
  gc()
  if(!is.null(cluster)) {
  # now, what about multiway clustering?
  # we can't scale with the residuals, but with the cluster demeaned residuals
  # and the group means
      skel <- lapply(cluster,function(cl) rep(0,nlevels(cl)))
      cdres <- demeanlist(smpdraw,cluster)
      rawmeans <- relist(kaczmarz(cluster,smpdraw-cdres),skel)
      cmeans <- lapply(rawmeans, function(cm) cm - mean(cm))
  } 
  for(i in 1:blks) {
    if(robust) {
      # robust residuals, variance is each squared residual
      rsamp <- rnorm(vpb*length(smpdraw))*abs(smpdraw)
    } else if(!is.null(cluster)) {
        # We draw group specific errors, sum them, and add a draw of within group specific errors
        rsamp <- as.vector(replicate(vpb, {
            rowSums(as.matrix(sapply(seq_along(cluster), function(cl) {
                clu <- cluster[[cl]]
                (cmeans[[cl]]*rnorm(nlevels(clu)))[clu]
            }))) + cdres*rnorm(length(R))
        }))
    } else {
      # IID residuals
      rsamp <- sample(smpdraw,vpb*length(smpdraw),replace=TRUE)
    }
    dim(rsamp) <- c(length(smpdraw),vpb)
    v <- kaczmarz(obj$fe,demeanlist(rsamp,obj$fe,eps=eps,threads=threads,means=TRUE)+Rvec, eps, threads=threads)*sefact
#    newR <- rsamp - demeanlist(rsamp,obj$fe,eps=eps,threads=threads) + Rvec
#    v <- kaczmarz(obj$fe,newR,eps,threads=threads)*sefact
#    rm(newR)
    rm(rsamp)
    efv <- apply(v,2,ef,addnames=FALSE)
    vsum <- vsum + rowSums(efv)
    vsq <- vsq + rowSums(efv**2)
    now <- as.integer(Sys.time())
    if(now-last > 300) {
      cat('...finished',i*vpb,'of',newN,'vectors in',now-start,'seconds\n')
      last <- now
    }
  }
  if(robust)
      sename <- 'robustse'
  else if(!is.null(cluster))
      sename <- 'clusterse'
  else
      sename <- 'se'
  if(!is.null(lhs)) sename <- paste(sename,lhs,sep='.')
  alpha[,sename] <- sqrt(vsq/newN - (vsum/newN)**2)/(1-0.75/newN-7/32/newN**2-9/128/newN**3)

  return(alpha)
}
