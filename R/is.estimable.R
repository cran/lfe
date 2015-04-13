# $Id: is.estimable.R 1655 2015-03-18 18:51:06Z sgaure $
# check whether a function is estimable
is.estimable <- function(ef,fe,R=NULL,nowarn=FALSE,keepdiff=FALSE, threshold=1e-5) {
  if(!is.function(ef)) stop('ef must be a function')
  
  N <- sum(unlist(lapply(fe,function(f) {
    x <- attr(f,'x')
    if(is.matrix(x)) nlevels(f)*ncol(x) else nlevels(f)
  })))
  if(is.null(R)) {
    # make a suitable residual
    nr <- length(fe[[1]])
    vec <- unlist(lapply(fe,function(f) {
      x <- attr(f,'x')
      if(is.matrix(x)) return(unlist(apply(x,2,function(cl) cl*runif(nlevels(f))[f])))
      r <- runif(nlevels(f))[f]
      if(is.null(x)) r else unlist(r*x)
    }))
    dim(vec) <- c(nr, length(vec)/nr)
    R <- rowSums(vec)
  }
  v1 <- ef(kaczmarz(fe,R,init=runif(N)),TRUE)
  v2 <- ef(kaczmarz(fe,R,init=runif(N)),TRUE)
  df <- max(abs(v1-v2))
  if(df > threshold) {
    bad <- which.max(abs(v1-v2))
    badname <- names(bad)
    if(!nowarn)
      warning('non-estimable function, largest error ',
              format(df,digits=1),' in coordinate ',bad, ' ("',badname,'")')
    return(structure(FALSE,diff=if(!keepdiff) NULL else v1-v2))
  }
  structure(TRUE,diff=if(!keepdiff) NULL else v1-v2)
}
