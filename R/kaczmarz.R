
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
  mat <- t(do.call(rBind,lapply(fl,as,'sparseMatrix')))
  if(is.null(dim(R))) dim(R) <- c(length(R),1)
# target function
  fn <- function(p,R) sum((mat %*% p - R)^2)
# its gradient
  gr <- function(p,R) {
      a <- 2*((mat %*% p) - R)
      dim(a) <- c(1,nrow(mat))
      as.vector(a %*% mat)
  }
  if(is.null(init)) init <- rnorm(ncol(mat))
  apply(R,2,function(vec) {
    cg <- Rcgmin::Rcgmin(init,fn,gr,control=list(eps=eps),R=vec)
    if(cg$convergence != 0) warning(cg$message)
    cg$par
  })
}

getfe.kaczmarz <- function(obj,se=FALSE,eps=getOption('lfe.eps'),ef='ref',bN=100,
                           robust=FALSE, cluster=NULL, lhs=NULL) {

  if(is.character(ef)) {
    ef <- efactory(obj,opt=ef)
  }
  if(!isTRUE(attr(ef,'verified')) && !is.estimable(ef, obj$fe)) {
    warning('Supplied function seems non-estimable')
  }

  multlhs <- length(obj$lhs) > 1
  if(is.null(lhs)) {
    R <- obj$r.residuals-obj$residuals
  } else {
    if(!all(lhs %in% obj$lhs)) stop('lhs must be subset of ', paste(obj$lhs, collapse=' '))
    R <- obj$r.residuals[,lhs, drop=FALSE] - obj$residuals[,lhs, drop=FALSE]
  }
  v <- kaczmarz(obj$fe,R,eps)

  if(is.matrix(v) && ncol(v) > 1) {
    v <- apply(v,2,ef,addnames=TRUE)
    vtmp <- ef(v[,1],addnames=TRUE)
    extra <- attr(vtmp, 'extra')
    nm <- names(vtmp)
  } else {
    v <- ef(v,TRUE)
    extra <- attr(v,'extra')
    nm <- names(v)
  }
  res <- data.frame(effect=v)
  if(multlhs) colnames(res) <- paste('effect',colnames(R),sep='.')
  if(!is.null(extra)) res <- cbind(res,extra)
  rownames(res) <- nm
  attr(res,'ef') <- ef

  if(se) {
    if(multlhs) {
      for(lh in colnames(R)) {
        res <- btrap(res,obj,bN,eps=eps, robust=robust, cluster=cluster, lhs=lh)        
      }
    } else {
      res <- btrap(res,obj,bN,eps=eps, robust=robust, cluster=cluster)
    }
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
# No, can't do that. The same level may occur in separate WW-components.
# WW is about partitioning of the dataset, not of the levels.

# return level names in appropriate order
# the G(x:f) with x a matrix makes it slightly complicated
xlevels <- function(n,f,sep='.') {
  x <- attr(f,'x')
  plev <- paste(n,levels(f),sep=sep)
  if(is.null(x) || !is.matrix(x)) return(plev)
  nam <- attr(f,'xnam')
  if(!is.matrix(x)) return(paste(nam,plev,sep=sep))
  matnam <- colnames(x)
  if(is.null(matnam)) matnam <- paste(nam,1:ncol(x),sep='') else matnam <- paste(nam,matnam,sep='')
  plev <- sub('.*:','',plev)
  return(as.vector(t(outer(matnam,plev,paste,sep=':'))))
}

nxlevels <- function(n,f) {
  x <- attr(f,'x')
  plev <- rep(n,nlevels(f))
  if(is.null(x) || !is.matrix(x)) return(plev)
  nam <- attr(f,'xnam')
  matnam <- colnames(x)
  if(is.null(matnam)) matnam <- paste(nam,1:ncol(x),sep='') else matnam <- paste(nam,matnam,sep='')
  plev <- sub('.*:','',plev)
  return(as.vector(t(outer(matnam,plev,paste,sep=':'))))
}




