# $Id: waldtest.R 1655 2015-03-18 18:51:06Z sgaure $
waldtest <- function(object, R, r, type=c('default','iid','robust','cluster'), lhs=NULL, df1, df2) {
  if(inherits(object,'felm') && object$nostats) stop('No Wald test for objects created with felm(nostats=TRUE)')
  # We make a chi^2 to test whether the equation R theta = r holds.
  # The chi^2 is computed according to Wooldridge (5.34, 10.59).
  # I.e. a Wald test W = N*(beta' (R V^{-1} R')^{-1} beta) where beta = R theta - r
  # W is chi2 with length(r) df,
  # and V is th covariance matrix.

  # First, find V. It's in either object$vcv, object$robustvcv or object$clustervcv
  if(is.null(lhs) && length(object$lhs) > 1) {
    stop('Please specify lhs=[one of ',paste(object$lhs, collapse=','),']')
  }
  if(!is.null(lhs) && is.na(match(lhs, object$lhs)))
      stop('Please specify lhs=[one of ',paste(object$lhs, collapse=','),']')

  type <- type[1]
  if(identical(type,'default')) {
    if(is.null(object$clustervar))
        V <- vcov(object, type='iid', lhs=lhs)
    else
        V <- vcov(object, type='cluster', lhs=lhs)
  } else
      V <- vcov(object, type=type, lhs=lhs)

#  if(is.null(lhs) && length(object$lhs) == 1) lhs <- object$lhs
  cf <- coef(object)
  if(is.matrix(cf))
      nmc <- rownames(cf)
  else
      nmc <- names(cf)

  if(inherits(R,'formula') || is.call(R) || is.name(R)) {
    Rr <- formtoR(R, nmc)
    R <- Rr[,-ncol(Rr), drop=FALSE]
    r <- Rr[,ncol(Rr)]
  } else  if(!is.matrix(R)) {
    # it's not a matrix, so it's a list of parameters, either
    # names, logicals or indices
    if(is.null(R)) R <- nmc
    if(is.character(R)) {
      ev <- match('endovars', R)
      if(!is.na(ev)) {
        # replace with endogenous variables
        R <- c(R[-ev],object$endovars)
      }
      # did user specify any of the endogenous variables?
      fitvars <- paste('`',R,'(fit)`',sep='')
      fitpos <- match(fitvars,nmc)
      # replace those which are not NA
      noNA <- which(!is.na(fitpos))
      R[noNA] <- fitvars[noNA]
      Ri <- match(R, nmc)
      if(anyNA(Ri)) stop("Couldn't find variables ",paste(R[is.na(Ri)],collapse=','))
      R <- Ri
    } else if(is.logical(R)) {
      R <- which(R)
    }
    # here R is a list of positions of coefficients
    # make the projection matrix.

    RR <- matrix(0,length(R),length(coef(object,lhs=lhs)))
    for(i in seq_along(R)) {
      RR[i,R[i]] <- 1
    }
    R <- RR
  } 
  if(missing(r) || is.null(r))
      r <- rep(0,nrow(R))
  else if(length(r) != nrow(R)) stop('nrow(R) != length(r)')
  cf <- coef(object, lhs=lhs)
  cf[is.na(cf)] <- 0
  beta <- R %*% cf - r
  V[is.na(V)] <- 0   # ignore NAs
  W <- try(sum(beta * solve(R %*% V %*% t(R),beta)))

  if(inherits(W,'try-error')) 
      W <- as.numeric(t(beta) %*% pinvx(R %*% V %*% t(R)) %*% beta)

  # W follows a chi2(Q) distribution, but the F-test has another
  # df which is ordinarily object$df. However, if there are clusters
  # the df should be reduced to the number of clusters-1
  if(missing(df2)) {
    df2 <- object$df
    if((!is.null(object$clustervar) && type %in% c('default','cluster')) ) {
      df2 <- min(nlevels(object$clustervar[[1]])-1, df2)
    }
  }

  if(missing(df1))
      df1 <- length(beta)

  F <- W/df1
  # F follows a F(df1,df2) distribution

  structure(c(p=pchisq(W, df1, lower.tail=FALSE), chi2=W, df1=df1,
              p.F=pf(F,df1,df2, lower.tail=FALSE), F=F, df2=df2),
            formula=Rtoform(R,r,nmc))
}

# convert a formula which is a set of linear combinations like ~x+x3 | x2-x4+3 to
# matrices R and r such that R %*% coefs = r
# the vector r is return as the last column of the result
formtoR <- function(formula, coefs) {

  conv <- function(f) formtoR(f, coefs)

  lf <- as.list(formula)
  if(lf[[1]] == as.name('~') || lf[[1]] == as.name('quote')) return(conv(lf[[2]]))
  # here we have a single formula w/o '~' in front, e.g. x+x3|x2-x4, or just x+x3
  # split off parts '|' in a loop
  R <- NULL
#  if(length(lf) != 1) stop('length of ',lf, ' is != 1')
#  lf <- as.list(lf[[1]])
  op <- lf[[1]]
  if(op == as.name('|')) {
    return(rbind(conv(lf[[2]]), conv(lf[[3]])))
  } else if(op == as.name('+')) {
    if(length(lf) == 2) return(conv(lf[[2]])) # unary +
    return(conv(lf[[2]]) + conv(lf[[3]]))
  } else if(op == as.name('-')) {
    if(length(lf) == 2) return(-conv(lf[[2]])) # unary -
    return(conv(lf[[2]]) - conv(lf[[3]]))    
  } else if(op == as.name('*')) {
    f1 <- conv(lf[[2]])
    f2 <- conv(lf[[3]])
    # the first one must be a numeric, i.e. only last column filled in
    # and it's negative
    fac <- -f1[length(f1)]
    return(fac * conv(lf[[3]])) 
  } else if(is.name(op)) {
    res <- matrix(0,1,length(coefs)+1)
    pos <- match(as.character(op), coefs)
    if(is.na(pos)) {
      ivspec <- paste("`",as.character(op),"(fit)`", sep='')
      pos <- match(ivspec, coefs)
    }
    if(is.na(pos)) stop("Can't find ", op, " among coefficients ", paste(coefs, collapse=','))
    res[pos] <- 1
    return(res)
  } else if(is.numeric(op)) {
    return(matrix(c(rep(0,length(coefs)), -op), 1))
  } else {
    stop('Unkwnown item ',as.character(op), ' in formula ',formula)
  }
}

Rtoform <- function(R,r, coefs) {
  coefs <- gsub('`','',coefs,fixed=TRUE)
  form <- paste('~',paste(apply(R, 1, function(row) {
    w <- which(row != 0)
    rw <- paste(' ', row[w], '*`', coefs[w], '`', collapse=' + ', sep='')
    rw <- gsub('+ -',' - ',rw, fixed=TRUE)
    rw <- gsub(' 1*','',rw, fixed=TRUE)
    rw <- gsub('(fit)','',rw, fixed=TRUE)
    rw
  }), ' + ', -r, collapse='|', sep=''))
  form <- gsub('+ -','-',form, fixed=TRUE)
  form <- gsub(' 0.',' .',form, fixed=TRUE)
  form <- gsub('+ 0','',form, fixed=TRUE)
  local(as.formula(form))
}
