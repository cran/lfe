# $Id: compfactor.R 1655 2015-03-18 18:51:06Z sgaure $
compfactor <- function(fl, WW=FALSE) {
  if(length(fl) == 0) return(factor(NULL))
  N <- length(fl[[1]])
  purefls <- sapply(fl,function(f) is.null(attr(f,'x')))
  fl <- fl[purefls]
  if(length(fl) <= 1) return(factor(rep(1,N)))
  if(WW && length(fl) > 2) {
    cf <- factor(.Call(C_wwcomp,fl))
  } else {
    cf <- factor(.Call(C_conncomp,fl[1:2]))
  }
  cf
}
