# Author: Simen Gaure
# Copyright: 2011, Simen Gaure
# Licence: Artistic 2.0
# pseudo-inverse, snatched from MASS

# this call changes or adds the dimnames of obj without duplicating it
# it should only be used on objects with a single reference
# Our use of it is safe, and we don't export it.
setdimnames <- function(obj, nm) {
  .Call(C_setdimnames,obj,nm)
}

pinvx <- function(X) {
#    return(pinv(nazero(X)))
    ch <- cholx(nazero(X))
    badvars <- attr(ch,'badvars')
    inv1 <- chol2inv(ch)
    if(is.null(badvars)) return(inv1)
    inv <- matrix(0,nrow(X),ncol(X))
    inv[-badvars,-badvars] <- inv1
    structure(inv,badvars=attr(ch,'badvars'))
}

# Do a Cholesky to detect multi-collinearities
cholx <- function(mat, eps=1e-6) {
  if(is.null(dim(mat))) dim(mat) <- c(1,1)
  N <- dim(mat)[1]
  if(N == 1) {
      return(structure(sqrt(mat),badvars=if(mat<=0) 1 else NULL))
  }

  # first, try a pivoted one
  tol <- N*getOption('lfe.eps')
#  tol <- -1
  chp <- chol(mat,pivot=TRUE,tol=tol)
  rank <- attr(chp,'rank')
  if(rank == N) return(chol(mat))
  pivot <- attr(chp,'pivot')
  oo <- order(pivot)
  badvars <- pivot[((rank+1):N)]
  ok <- (1:N)[-badvars]
  ch <- chol(mat[ok,ok])
  return(structure(ch,badvars=badvars))
}


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

nazero <- function(x) ifelse(is.na(x),0,x)

ftest <- function(z, zr=NULL, vcov=z$vcv) {
  rdf <- z$df #z$N - z$p - 1
  if(is.null(zr)) {
# we should do F-test vs. model with intercept.
# but the intercept in z is implicit.
    F <- (t(coef(z)) %*% pinvx(vcov) %*% coef(z))/z$p
    return( c(F=F, p=pf(F, z$p, rdf, lower.tail=FALSE), df1=z$p, df2=rdf))
  }
#  df1 <- length(coef(z)) - length(coef(zr))
  df1 <- z$p - zr$p
  c1 <- coef(z)
  c2 <- rep(0,length(c1))
  names(c2) <- names(c1)
  c2[names(coef(zr))] <- coef(zr)
  F <- (t(c1-c2) %*% pinvx(vcov) %*% (c1-c2))/df1
  return(c(F=F, p=pf(F,df1, rdf, lower.tail=FALSE), df1=df1,df2=rdf))
}


mkgraph <- function(f1,f2) {
  if(!requireNamespace('igraph', quietly=TRUE)) stop('Package igraph not found')
#  graph.edgelist(cbind(paste('f1',f1),paste('f2',f2)), directed=FALSE)
#  graph.edgelist(cbind(500000000+as.integer(f1),f2), directed=FALSE)
  igraph::graph.adjacency(tcrossprod(do.call(rBind,
                                     lapply(list(f1,f2), as, 'sparseMatrix')))>0,
                  'undirected', diag=FALSE)
}

diamgraph <- function(flist,approx=TRUE) {
  if(!requireNamespace('igraph', quietly=TRUE)) stop('Package igraph not found')
  gr <- mkgraph(flist[[1]],flist[[2]])
# find largest cluster
  cl <- igraph::clusters(gr)$membership
  lcl <- which(cl == which.max(table(cl)))
  if(approx) 
    max(igraph::shortest.paths(gr,v=sample(lcl,10),to=sample(lcl,10)))
  else
    igraph::diameter(igraph::induced.subgraph(gr,lcl))
}

diammatrix <- function(flist, approx=TRUE) {

  flen <- length(flist)
  if(flen < 2) return(0)
  val <- matrix(0,flen,flen)
  colnames(val) <- names(flist)
  rownames(val) <- names(flist)
  for(if1 in 1:(flen-1))
    for(if2 in (if1+1):flen) {
      val[if1,if2] <- diamgraph(flist[c(if1,if2)], approx)
    }

# fill in lower triangle:
  val[row(val) > col(val)] <- t(val)[row(val) > col(val)]
  val
}

# if(!exists('fixef')) fixef <- function(object,...) UseMethod('fixef')

# compute rank deficiency of D-matrix
rankDefic <- function(fl,method='cholesky',mctol=1e-3) {
  eps <- sqrt(.Machine$double.eps)
  if(length(fl) == 1) return(1)
  if(length(fl) == 2) return(nlevels(compfactor(fl)))
  if(method == 'cholesky') {
    D <- makeDmatrix(fl)
    Ch <- as(Cholesky(crossprod(D), super=TRUE, perm=TRUE, Imult=eps), 'sparseMatrix')
    sum(diag(Ch) < eps^(1/4))
  } else if(method == 'mc') {
    totlev <- sum(sapply(fl,nlevels))
    N <- length(fl[[1]])
    len <- length(fl) + 1
    # now, we will use the rank deficiency d to compute
    # degrees of freedom corrections. I.e. (tr+totlev)+len should
    # be within a certain relative tolerance, which translates
    # to an absolute tolerance of tr
    tolfun <- function(tr) -mctol*abs((tr+totlev)+len)
    initr <- 0
    N - (mctrace(fl,tol=tolfun,initr=initr) + totlev)+len
    
  } else {
    D <- makeDmatrix(fl)
    as.integer(ncol(D) - rankMatrix(crossprod(D), method='qr.R'))
  }
}


makeDmatrix <- function(fl) {
  # make the D matrix
  # for pure factors f, it's just the t(as(f,'sparseMatrix'))
  # if there's a covariate vector x, it's t(as(f,'sparseMatrix'))*x
  # for a covariate matrix x, it's the cbinds of the columns with the factor-matrix
  do.call(cBind,lapply(fl, function(f) {
    x <- attr(f,'x')
    fm <- t(as(f,'sparseMatrix'))
    if(is.null(x)) return(fm)
    if(!is.matrix(x)) return(fm*x)
    do.call(cBind,apply(x,2,'*',fm))
  }))
}

makePmatrix <- function(fl) {
  D <- makeDmatrix(fl)
  DtD <- crossprod(D)
  DtDi <- pinvx(as.matrix(DtD))
  badvars <- attr(DtDi,'badvars')
  if(is.null(badvars)) return(D %*% DtDi %*% t(D))
  D <- D[,-badvars,drop=FALSE]
  return(D %*% pinvx(as.matrix(crossprod(D))) %*% t(D))
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
  if(missing(cf)) cf <- compfactor(fl)
  numpure <- sum(sapply(fl,function(f) is.null(attr(f,'x'))))
  if(numpure == 1) return(0)
  if(numpure == 2) return(nlevels(cf))
  if(identical(exactDOF,'rM')) {
    return(rankDefic(fl, method='qr'))
  } else if(identical(exactDOF,'mc')) {
    return(rankDefic(fl, method='mc'))
  } else  if(exactDOF) {
    return(rankDefic(fl, method='cholesky'))
  }
  return(nlevels(cf) + numpure-2)
}

#prettyprint a list of integers
ilpretty <- function(il) {
  paste(tapply(il,cumsum(c(1,diff(il)) != 1),
         function(x) if(length(x) == 1) {
           as.character(x)
         } else {
           paste(x[[1]],x[[length(x)]],sep=':')
         }), collapse=' ')
}

makefitnames <- function(s) paste('`',s,'(fit)`',sep='')
delete.icpt <- function(x) {
  asgn <- attr(x,'assign')
  icpt <- asgn == 0
  if(!length(icpt)) return(x)
  asgn <- asgn[!icpt]
  ctr <- attr(x,'contrasts')
  x <- x[,!icpt,drop=FALSE]
  attr(x,'assign') <- asgn
  attr(x,'contrasts') <- ctr
  x
}
