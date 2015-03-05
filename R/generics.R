print.felm <- function(x,digits=max(3,getOption('digits')-3),...) {
  z <- x
  if(z$p == 0) {
    cat('(No coefficients)\n')
    return()
  }
  print(coef(x),digits=digits,...)
}

#fixef.felm <- function(object,...) {
#  fe <- getfe(object,...)
#  f <- fe[,'fe']
#  l <- lapply(levels(f),function(n) {v <- fe[f == n,'effect']; names(v) <- as.character(fe[f==n,'idx']); v})
#  names(l) <- levels(f)
#  l
#}

coef.felm <- function(object, ..., lhs=NULL) {
  if(is.null(lhs)) {
    if(ncol(object$coefficients) == 1) return({r <- as.vector(object$coefficients); names(r) <- rownames(object$coefficients); r})
    object$coefficients
  } else {
#    if(anyNA(match(lhs, colnames(object$coefficients))))
    if(any(!(lhs %in% colnames(object$coefficients))))
        stop('Please specify lhs as one of ',paste(object$lhs, collapse=','))
    object$coefficients[,lhs,drop=FALSE]
  }
}

residuals.felm <- function(object, ..., lhs=NULL) {
  if(is.null(lhs)) {
    if(ncol(object$coefficients) == 1) return({r <- as.vector(object$coefficients); names(r) <- rownames(object$coefficients); r})
    object$coefficients
  } else {
#    if(anyNA(match(lhs, colnames(object$coefficients))))
    if(any(!(lhs %in% colnames(object$coefficients))))
        stop('Please specify lhs as one of ',paste(object$lhs, collapse=','))
    object$residuals[,lhs,drop=FALSE]
  }
}

vcov.felm <- function(object,...,type=c('iid','robust','cluster'),lhs=NULL) {
#  if(is.na(match(type[1], c('iid', 'robust', 'cluster'))))
  if(!(type[1] %in% c('iid', 'robust', 'cluster')))
      stop("specify vcov-type as 'iid', 'robust' or 'cluster'")

  if(is.null(lhs) && length(object$lhs) > 1) {
    stop('Please specify which lhs to retrieve vcov for with vcov(...,lhs=[one of ',
         paste(object$lhs,collapse=','),'])')
  }
  if(is.null(lhs) || length(object$lhs) == 1) {
    if(type[1] == 'iid') return(object$vcv)
    if(type[1] == 'robust') return(object$robustvcv)
    if(type[1] == 'cluster') return(object$clustervcv)
  }

  if(is.na(match(lhs, object$lhs))) {
    stop('Please specify which lhs to retrieve vcov for with vcov(...,lhs=[one of ',
         paste(object$lhs,collapse=','),'])')
  }
  if(type[1] == 'iid') return(object$STATS[[lhs]]$vcv)
  if(type[1] == 'robust') return(object$STATS[[lhs]]$robustvcv)
  if(type[1] == 'cluster') return(object$STATS[[lhs]]$clustervcv)
}

update.felm <- function (object, formula., ..., evaluate = TRUE) 
{
    if (is.null(call <- getCall(object))) 
        stop("need an object with call component")

    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.)) 
        call$formula <- formula(update(as.Formula(call$formula), formula.))
    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate) 
        eval(call, parent.frame())
    else call
}

xtable.summary.felm <- function(x, caption=NULL, label=NULL, align=NULL, digits=NULL,
                                display=NULL, ...) {
    cl <- match.call(expand.dots=FALSE)
    do.call(getS3method('xtable','summary.lm'), as.list(cl)[-1])
}


xtable.felm <- function(x, caption=NULL, label=NULL, align=NULL, digits=NULL,
                        display=NULL, ...) {
    cl <- match.call(expand.dots=FALSE)
    do.call(getS3method('xtable','lm'), as.list(cl)[-1])
}

print.summary.felm <- function(x,digits=max(3L,getOption('digits')-3L),...) {
  if(!is.null(x$lhs))
      cat('Summary for outcome',x$lhs,'\n')
  cat('\nCall:\n  ',deparse(x$call),'\n\n')

  qres <- zapsmall(quantile(x$residuals), digits+1L)
  cat('Residuals:\n')
  names(qres) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(qres,digits=digits,...)

  cat('\nCoefficients:\n')
  if(x$Pp <= 0)
    cat("(No coefficients)\n")
  else {
    printCoefmat(x$coefficients,digits=digits)
    cat('\nResidual standard error:',format(signif(x$rse,digits)),'on',x$rdf,'degrees of freedom\n')
    if (nzchar(mess <- naprint(x$na.action)))
        cat("  (", mess, ")\n", sep = "")
    cat('Multiple R-squared(full model):',formatC(x$r2,digits=digits),'  Adjusted R-squared:',
        formatC(x$r2adj,digits=digits),'\n')
    if(!is.null(x$P.r.squared)) {
      cat('Multiple R-squared(proj model):',formatC(x$P.r.squared,digits=digits),
          '  Adjusted R-squared:',  formatC(x$P.adj.r.squared,digits=digits),'\n')
    }
    if(x$badF)
      cat('F-statistic(full model, *iid*):')
    else
      cat('F-statistic(full model):')
    hasicpt <- if(is.null(x$hasicpt)) 0 else 1

    if(is.null(x$F.fstat)) {
      cat(formatC(x$fstat,digits=digits),'on',x$p-hasicpt,'and',x$rdf,
          'DF, p-value:',format.pval(x$pval,digits=digits),'\n')
    } else {
      cat(formatC(x$F.fstat['F'],digits=digits),'on',x$F.fstat['df1'],
          'and',x$F.fstat['df2'], 'DF, p-value:',
          format.pval(x$F.fstat['p'],digits=digits),'\n')
    }
    cat('F-statistic(proj model):',formatC(x$P.fstat[['F']],digits=digits),'on',
        x$P.fstat[['df1']], 'and',x$P.fstat[['df2']], 'DF, p-value:',
        format.pval(x$P.fstat[['p.F']],digits=digits),'\n')

    if(!is.null(x$iv1fstat)) {
      if1 <- x$iv1fstat
      cat('F-statistic(excl instr.):')
      
      cat(formatC(if1[['F']],digits=digits),'on',
          if1[['df1']],'and',if1[['df2']],'DF, p-value:',
          format.pval(if1[['p.F']],digits=digits),'\n')      
    }
    
    if(!is.null(x$E.fstat)) {
      if1 <- x$E.fstat
      cat('F-statistic(endog. vars):')
      cat(formatC(if1[['F']],digits=digits),'on',
          if1[['df1']],'and',if1[['df2']],'DF, p-value:',
          format.pval(if1[['p.F']],digits=digits),'\n')      

    }

    if(length(x$fe) > 2 && !identical(x$exactDOF,'rM') && !x$exactDOF)
        cat('*** Standard errors may be too high due to more than 2 groups and exactDOF=FALSE\n')
  }
  if(!is.null(x$numctrl) && x$numctrl != 0)
      cat(x$numctrl, 'variable(s) were projected out\n')
  cat('\n\n')
}

summary.felm <- function(object,...,robust=!is.null(object$clustervar),lhs=NULL) {
  z <- object
  if(z$nostats) stop('No summary for objects created with felm(nostats=TRUE)')
  if(is.null(lhs)) {
    if(length(z$lhs) > 1)
        stop('Please specify lhs=[one of ',paste(z$lhs,collapse=','),']')
    STATS <- z
    lhs <- object$lhs
  } else {
    if(is.na(match(lhs, z$lhs))) 
        stop('Please specify lhs=[one of ',paste(z$lhs,collapse=','),']')
    if(length(z$lhs) >= 1)
        STATS <- z$STATS[[lhs]]
    else
        STATS <- z
  }
  res <- list()
  res$p <- z$p
  res$Pp <- z$Pp
  res$numctrl <- z$numctrl
  res$ctrlnames <- z$ctrlnames
  if(length(z$lhs) > 1) res$lhs <- lhs
  if(res$Pp == 0) {
    res <- list(residuals=as.vector(z$residuals[,lhs]),p=z$p,Pp=0,call=z$call)
    class(res) <- 'summary.felm'
    return(res)
  }

  res$terms <- z$terms
  res$call <- z$call
  res$badF <- FALSE
  if(is.logical(robust) && robust) {
    if(!is.null(STATS$cse)) {
      coefficients <- cbind(z$beta[,lhs],STATS$cse,STATS$ctval,STATS$cpval)
      sdnam <- 'Cluster s.e.'
      res$badF <- TRUE
    } else {
      coefficients <- cbind(z$beta[,lhs],STATS$rse,STATS$rtval,STATS$rpval)
      sdnam <- 'Robust s.e'
      res$badF <- TRUE
    }
  } else {
    sdnam <- 'Std. Error'
    coefficients <- cbind(z$beta[,lhs],STATS$se,STATS$tval,STATS$pval)
  }

  if(!is.null(coefficients)) {
    dimnames(coefficients) <- 
       list(rownames(z$beta),
           c('Estimate',sdnam,'t value','Pr(>|t|)'))

  }
  res$coefficients <- coefficients
  res$residuals <- as.vector(z$residuals[,lhs])

  qres <- quantile(res$residuals,na.rm=TRUE)
  names(qres) <- c("Min", "1Q", "Median", "3Q", "Max")
  res$qres <- qres

  # compute
  # residual se, df
  # mrsq, adj rsq
  # fstat, p-value

  p <- z$p
  if(is.null(z$hasicpt)) hasicpt <- 0 else hasicpt <- z$hasicpt
  if(length(z$fe) > 0) hasicpt <- 1
  res$hasicpt <- hasicpt

  rdf <- z$N - p
  rss <- sum(z$residuals[,lhs]^2)

  if(!is.null(z$TSS))
      tss <- z$TSS[[lhs]]
  else
      tss <- sum( (z$response[,lhs] - mean(z$response[,lhs]))^2)

  mss <- tss - rss

#  mss2 <- if(hasicpt) sum((z$fitted.values[,lhs]-mean(z$fitted.values[,lhs]))^2) else sum(z$fitted.values[,lhs]^2)
#  tss <- mss + rss

  resvar <- rss/rdf
  r2 <- mss/tss
  r2adj <- 1-(1-r2)*((z$N-hasicpt)/rdf)
  
# F-tests for 2. stage iv is different
  if(!is.null(z$iv.residuals)) {
    # We have F = (tss - rss)/rss  (and some df factor)
    # however, the numerator should be residuals w.r.t. to the
    # fitted variables whereas the denominator should be w.r.t. to
    # the original variables. (Wooldridge, p. 99)
    # every metric verified to match Stata ivregress with small sample adjustment Jan 12, 2015
    mss <- tss - sum(z$iv.residuals[,lhs]^2)
  }
  # hmm, fstat should be computed differently when s.e. are robust or clustered.

# full model Fstat for i.i.d
  F <- as.numeric((mss/(p-hasicpt))/resvar)  # get rid of name
  res$F.fstat <- c(F=F,
                   df1=p-hasicpt,
                   df2=rdf,
                   p=pf(F,p-hasicpt,rdf,lower.tail=FALSE))

  res$fstat <- F
  res$pval <- res$F.fstat['p']
  # projected R2?  I.e. how much is explained by the variables that
  # have not been projected out?  That's similar to the projected F-test.
  # So the tss is not the tss of the response variable, but of the
  # projected response  (P.response doesn't exist in the old felm)
  if(!is.null(z$P.TSS)) {
    tss <- z$P.TSS[lhs]
    mss <- tss - rss
    res$P.r.squared <- mss/tss
    res$P.adj.r.squared <- 1-(1-res$P.r.squared)*((z$N-hasicpt)/rdf)
  }
  
  # use wald test if robust or clustered
  mcoef <- rownames(z$coefficients)
  mcoef <- mcoef[!(mcoef %in% '(Intercept)')]

  wtype <- 'iid'
  if(robust) wtype <- if(is.null(z$clustervar)) 'robust' else 'cluster'
  F <- try(waldtest(z,mcoef, type=wtype, lhs=lhs))
  if(inherits(F,'try-error')) {
    warning("can't compute cluster F-test")
    F <- list(F=NaN,p.F=NaN)
  }

  res$P.fstat <- F

  # then a Wald test on the endogenous variables
  if(!is.null(z$iv.residuals)) {
    res$E.fstat <- waldtest(z, 'endovars', type=wtype, lhs=lhs)
  }

  sigma <- sqrt(resvar)  

  res$exactDOF <- z$exactDOF
  if(is.list(z$iv1fstat)) {
    res$iv1fstat <- z$iv1fstat[[lhs]]
    res$rob.iv1fstat <- z$rob.iv1fstat[[lhs]]
  } else {
    res$iv1fstat <- z$iv1fstat
    res$rob.iv1fstat <- z$rob.iv1fstat
  }
  res$df <- c(rdf,rdf)
  res$sigma <- res$rse <- sigma
  res$rdf <- rdf
  res$r.squared <- res$r2 <- r2
  res$adj.r.squared <- res$r2adj <- r2adj
  res$fe <- z$fe
  res$N <- z$N
  res$na.action <- z$na.action
  class(res) <- 'summary.felm'
  res
}

