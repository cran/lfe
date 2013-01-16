# Author: Simen Gaure
# Copyright: 2011, Simen Gaure
# Licence: Artistic 2.0


summary.felm <- function(object,...,robust=FALSE) {
  z <- object
  res <- list()
  res$p <- z$p
  if(res$p == 0) {
    res <- list(residuals=as.vector(z$full.residuals),p=0,call=z$call)
    class(res) <- 'summary.felm'
    return(res)
  }
  res$terms <- z$terms
  res$call <- z$call
  if(is.logical(robust) && robust) {
    if(!is.null(object$cse)) {
      coefficients <- cbind(z$coefficients,z$cse,z$ctval,z$cpval)
      sdnam <- 'Cluster s.e.'
    } else {
      coefficients <- cbind(z$coefficients,z$rse,z$rtval,z$rpval)
      sdnam <- 'Robust s.e'
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
  res$residuals <- as.vector(z$full.residuals)

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
  rss <- sum(z$full.residuals^2)

  resvar <- rss/rdf
  sigma <- sqrt(resvar)

  tss <- sum( (z$response - mean(z$response))^2)
  mss <- tss - rss

  r2 <- mss/tss
  r2adj <- 1-(1-r2)*(z$N/rdf)
  fstat <- (mss/p)/resvar

  pval <- pf(fstat,p,rdf,lower.tail=FALSE)
  res$df <- c(rdf,rdf)
  res$rse <- sigma
  res$rdf <- rdf
  res$r2 <- r2
  res$r2adj <- r2adj
  res$fstat <- fstat
  res$pval <- pval
  res$fe <- z$fe
  class(res) <- 'summary.felm'
  res
}

print.summary.felm <- function(x,digits=max(3,getOption('digits')-3),...) {
#  if(x$p == 0) {
#    cat('(No coefficients)')
#    return()
#  }

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
    cat('Multiple R-squared:',formatC(x$r2,digits=digits),'  Adjusted R-squared:',
        formatC(x$r2adj,digits=digits),'\n')
    cat('F-statistic:',formatC(x$fstat,digits=digits),'on',x$p,'and',x$rdf,'DF, p-value:',format.pval(x$pval,digits=digits),'\n')
    if(length(x$fe) > 2)
      cat('*** Standard errors may be slightly too high due to more than 2 groups\n')
  }
  cat('\n\n')
  
}

print.felm <- function(x,digits=max(3,getOption('digits')-3),...) {
  z <- x
  if(z$p == 0) {
    cat('(No coefficients)\n')
    return()
  }
  print(coef(x),digits=digits)
  return()
  ans <- z[c('call','terms')]
  ans$residuals <- z$full.residuals
  rdf <- z$df
  ans$coefficients <- cbind(z$coefficients,z$se,z$tval,z$pval)
  dimnames(ans$coefficients) <- list(names(z$coefficients),
                                     c('Estimate','Std. Error','t value','Pr(>|t|)'))
  cat('\nCall:\n  ',deparse(z$call),'\n\n')
  cat('Residuals:\n')
  qres <- quantile(ans$residuals)
  names(qres) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(qres,digits=digits)

  cat('\nCoefficients:\n')
  printCoefmat(ans$coefficients,digits=digits)


  # compute
  # residual se, df
  # mrsq, adj rsq
  # fstat, p-value
  p <- z$p
  rdf <- z$N - p
  rss <- sum(z$full.residuals^2)
  f <- z$fitted
  resvar <- rss/rdf
  sigma <- sqrt(resvar)
#  mss <- sum( (f - mean(f))^2)
  mss <- sum(f^2)
  r2 <- mss/(mss+rss)
  r2adj <- 1-(1-r2)*(z$N/rdf)
  fstat <- (mss/p)/resvar
#  z <- list(call=z$call,residuals=z$full.residuals,
#            df=z$df)
  pval <- pf(fstat,p,rdf,lower.tail=FALSE)
  cat('\nResidual standard error:',format(signif(sigma,digits)),'on',rdf,'degrees of freedom\n')
  cat('Multiple R-squared:',formatC(r2,digits=digits),'  Adjusted R-squared:',
      formatC(r2adj,digits=digits),'\n')
  cat('F-statistic:',formatC(fstat,digits=digits),'on',p,'and',rdf,'DF, p-value:',format.pval(pval,digits=digits),'\n')
  if(length(z$fe) > 2)
    cat('*** Standard errors are slightly too high due to more than 2 groups\n')
  cat('\n\n')

}

fixef.felm <- function(object,...) {
  fe <- getfe(object,...)
  f <- fe[,'fe']
  l <- lapply(levels(f),function(n) {v <- fe[f == n,'effect']; names(v) <- as.character(fe[f==n,'idx']); v})
  names(l) <- levels(f)
  l
}

# if(!exists('fixef')) fixef <- function(object,...) UseMethod('fixef')
