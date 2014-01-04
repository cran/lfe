# Author: Simen Gaure
# Copyright: 2011, Simen Gaure
# Licence: Artistic 2.0
# pseudo-inverse, snatched from MASS


pinv <- function (X, tol = sqrt(.Machine$double.eps)) {
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X))
    X <- as.matrix(X)
  Xsvd <- svd(X)
  badvars <- integer(0)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  Positive <- ifelse(is.na(Positive),FALSE,Positive)
  badvars <- which(!Positive)
  if (all(Positive))
    res <- Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    res <- array(0, dim(X)[2L:1L])
  else
    res <- Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
                            t(Xsvd$u[, Positive, drop = FALSE]))
  structure(res,badvars=badvars)
}

ftest <- function(z, zr=NULL, vcov=z$vcv) {
  rdf <- z$df #z$N - z$p - 1
  if(is.null(zr)) {
# we should do F-test vs. model with intercept.
# but the intercept in z is implicit.
    F <- (t(coef(z)) %*% pinv(vcov) %*% coef(z))/z$p
    return( c(F=F, p=pf(F, z$p, rdf, lower.tail=FALSE), df1=z$p, df2=rdf))
  }
#  df1 <- length(coef(z)) - length(coef(zr))
  df1 <- z$p - zr$p
  c1 <- coef(z)
  c2 <- rep(0,length(c1))
  names(c2) <- names(c1)
  c2[names(coef(zr))] <- coef(zr)
  F <- (t(c1-c2) %*% pinv(vcov) %*% (c1-c2))/df1
  return(c(F=F, p=pf(F,df1, rdf, lower.tail=FALSE), df1=df1,df2=rdf))
}


summary.felm <- function(object,...,robust=!is.null(object$cse)) {
  z <- object
  res <- list()
  res$p <- z$p
  if(res$p == 0) {
    res <- list(residuals=as.vector(z$residuals),p=0,call=z$call)
    class(res) <- 'summary.felm'
    return(res)
  }
  res$terms <- z$terms
  res$call <- z$call
  res$badF <- FALSE
  if(is.logical(robust) && robust) {
    if(!is.null(object$cse)) {
      coefficients <- cbind(z$coefficients,z$cse,z$ctval,z$cpval)
      sdnam <- 'Cluster s.e.'
      res$badF <- TRUE
    } else {
      coefficients <- cbind(z$coefficients,z$rse,z$rtval,z$rpval)
      sdnam <- 'Robust s.e'
      res$badF <- TRUE
    }
  } else {
    sdnam <- 'Std. Error'
    coefficients <- cbind(z$coefficients,z$se,z$tval,z$pval)
  }
  if(!is.null(coefficients)) {
    dimnames(coefficients) <- 
       list(names(z$coefficients),
           c('Estimate',sdnam,'t value','Pr(>|t|)'))

  }
  res$coefficients <- coefficients
  res$residuals <- as.vector(z$residuals)

  qres <- quantile(res$residuals,na.rm=TRUE)
  names(qres) <- c("Min", "1Q", "Median", "3Q", "Max")
  res$qres <- qres

  # compute
  # residual se, df
  # mrsq, adj rsq
  # fstat, p-value

  p <- z$p

# should we subtract 1 for the intercept?
# a similar adjustment is done in summary.felm when computing p

  rdf <- z$N - p - 1
  rss <- sum(z$residuals^2)


  tss <- sum( (z$response - mean(z$response))^2)
  mss <- tss - rss

# things for 2. stage iv is different
  if(!is.null(z$ivresid)) {
    # residual from observed instrumented variables
    # residuals from null model - residuals from full model
    rs <- z$residuals - z$ivresid
    resvar <- sum(rs^2)/rdf
    r2 <- (tss-sum(rs^2))/tss
  } else {
    r2 <- mss/tss
    resvar <- rss/rdf
  }

  # hmm, fstat should be computed differently when s.e. are robust or clustered.
  

  res$fstat <- (mss/p)/resvar
  res$pval <- pf(res$fstat,p,rdf,lower.tail=FALSE)


  # use wald test if robust or clustered
  if(robust) {
    if(is.null(z$cse)) {
      F <- ftest(z,vcov=z$robustvcv)
    } else {
      F <- ftest(z,vcov=z$clustervcv)    
    }
    res$rob.fstat <- F[['F']]
    res$rob.pval  <- F[['p']]
  }

  sigma <- sqrt(resvar)  
  r2adj <- 1-(1-r2)*(z$N/rdf)


  res$exactDOF <- z$exactDOF
  if(!is.null(z$iv1fstat)) {
    res$iv1fstat <- z$iv1fstat
    res$iv1pval <- z$iv1fstat[['p']]
    if(robust) {
      if(is.null(z$cse))
        if1 <- z$rob.iv1fstat
      else
        if1 <- z$clu.iv1fstat
      res$rob.iv1fstat <- if1
      res$rob.iv1pval <- if1[['p']]
    }
  }

  res$df <- c(rdf,rdf)
  res$rse <- sigma
  res$rdf <- rdf
  res$r2 <- r2
  res$r2adj <- r2adj
#  res$fstat <- fstat
#  res$pval <- pval
  res$fe <- z$fe
  res$N <- z$N
  res$na.action <- z$na.action
  class(res) <- 'summary.felm'
  res
}

print.summary.felm <- function(x,digits=max(3,getOption('digits')-3),...) {
  cat('\nCall:\n  ',deparse(x$call),'\n\n')
  qres <- quantile(x$residuals)
  cat('Residuals:\n')
  names(qres) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(qres,digits=digits)

  cat('\nCoefficients:\n')
  if(x$p <= 0)
    cat("(No coefficients)\n")
  else {
    printCoefmat(x$coefficients,digits=digits)
    cat('\nResidual standard error:',format(signif(x$rse,digits)),'on',x$rdf,'degrees of freedom\n')
    if (nzchar(mess <- naprint(x$na.action)))
      cat("  (", mess, ")\n", sep = "")
    cat('Multiple R-squared:',formatC(x$r2,digits=digits),'  Adjusted R-squared:',
        formatC(x$r2adj,digits=digits),'\n')

    if(x$badF)
      cat('F-statistic(normal s.e.):')
    else
      cat('F-statistic:')
    cat(formatC(x$fstat,digits=digits),'on',x$p,'and',x$rdf,'DF, p-value:',format.pval(x$pval,digits=digits),'\n')
    if(x$badF) {
      cat('F-statistic(proj):',formatC(x$rob.fstat,digits=digits),'on',x$p,'and',x$rdf,'DF, p-value:',format.pval(x$rob.pval,digits=digits),'\n')
    }

    if(!is.null(x$iv1fstat)) {
      if1 <- x$iv1fstat
      if(x$badF)
        cat('F-stat (excl i.v. n.s.e):')
      else
        cat('F-stat (excl i.v.):')

      cat(formatC(if1[['F']],digits=digits),'on',
          if1[['df1']],'and',if1[['df2']],'DF, p-value:',
          format.pval(if1[['p']],digits=digits),'\n')      

      if(x$badF) {
        if1 <- x$rob.iv1fstat
        cat('F-stat (excl i.v. proj):',formatC(if1[['F']],digits=digits),'on',
          if1[['df1']],'and',if1[['df2']],'DF, p-value:',
          format.pval(if1[['p']],digits=digits),'\n')
      }
    }

    if(length(x$fe) > 2 && !identical(x$exactDOF,'rM') && !x$exactDOF)
      cat('*** Standard errors may be too high due to more than 2 groups and exactDOF=FALSE\n')
  }
  cat('\n\n')
  
}

vcov.felm <- function(object,...) return(object$vcv)

mkgraph <- function(f1,f2) {
  graph.edgelist(cbind(paste('f1',f1),paste('f2',f2)), directed=FALSE)
#  graph.edgelist(cbind(500000000+as.integer(f1),f2), directed=FALSE)
#  graph.adjacency(tcrossprod(do.call('rBind',
#                                     lapply(flist, as, 'sparseMatrix')))>0,
#                  'undirected', diag=FALSE)
}

diamgraph <- function(flist,approx=TRUE) {
  gr <- mkgraph(flist[[1]],flist[[2]])
# find largest cluster
  cl <- clusters(gr)$membership
  lcl <- which(cl == which.max(table(cl)))
  if(approx) 
    max(shortest.paths(gr,v=sample(lcl,10),to=sample(lcl,10)))
  else
    diameter(induced.subgraph(gr,lcl))
}

diammatrix <- function(flist, approx=TRUE) {
  flen <- length(flist)
  if(flen < 2) return(0)
  val <- matrix(0,flen,flen)
  colnames(val) <- names(flist)
  rownames(val) <- names(flist)
  for(if1 in 1:(flen-1))
    for(if2 in (if1+1):flen) {
      val[if1,if2] <- diamgraph(flist[c(if1,if2)])
    }

# fill in lower triangle:
  val[row(val) > col(val)] <- t(val)[row(val) > col(val)]
  val
}

print.felm <- function(x,digits=max(3,getOption('digits')-3),...) {
  z <- x
  if(z$p == 0) {
    cat('(No coefficients)\n')
    return()
  }
  print(coef(x),digits=digits)
}

fixef.felm <- function(object,...) {
  fe <- getfe(object,...)
  f <- fe[,'fe']
  l <- lapply(levels(f),function(n) {v <- fe[f == n,'effect']; names(v) <- as.character(fe[f==n,'idx']); v})
  names(l) <- levels(f)
  l
}

# if(!exists('fixef')) fixef <- function(object,...) UseMethod('fixef')

# compute rank deficiency of D-matrix
rankDefic <- function(fl,method='cholesky') {
  eps <- sqrt(.Machine$double.eps)
 # Ch <- as(Cholesky(tcrossprod(do.call('rBind',lapply(fl,as,Class='sparseMatrix'))),
#                  super=TRUE,perm=TRUE,Imult=eps),'sparseMatrix')
  D <- makeDmatrix(fl)
  if(method == 'cholesky') {
    Ch <- as(Cholesky(crossprod(D), super=TRUE, perm=TRUE, Imult=eps), 'sparseMatrix')
    sum(diag(Ch)^2 < eps^(1/3))
  } else {
    as.integer(ncol(D) - rankMatrix(D,method='qrLINPACK'))
  }
}


makeDmatrix <- function(fl) {
  # make the D matrix
  # for pure factors f, it's just the t(as(f,'sparseMatrix'))
  # if there's a covariate vector x, it's t(as(f,'sparseMatrix'))*x
  # for a covariate matrix x, it's the cbinds of the columns with the factor-matrix
  do.call('cBind',lapply(fl, function(f) {
    x <- attr(f,'x')
    fm <- t(as(f,'sparseMatrix'))
    if(is.null(x)) return(fm)
    if(!is.matrix(x)) return(fm*x)
    do.call('cBind',apply(x,2,'*',fm))
  }))
}

# total number of variables projected out
totalpvar <- function(fl) {
  if(length(fl) == 0) return(0)
  sum(sapply(fl, function(f) {
    x <- attr(f,'x')
    if(is.null(x) || !is.matrix(x)) return(nlevels(f))
    return(ncol(x)*nlevels(f))
  }))
}

nrefs <- function(fl, cf, exactDOF=FALSE) {
  if(length(fl) == 0) return(0)
  numpure <- sum(sapply(fl,function(f) is.null(attr(f,'x'))))
  if(identical(exactDOF,'rM')) {
    return(rankDefic(fl, method='qr'))
  } else  if(exactDOF) {
    return(rankDefic(fl, method='cholesky'))
  } else {
    return((if(numpure>1) nlevels(cf) else 0) + max(numpure-2, 0))
  }
}
