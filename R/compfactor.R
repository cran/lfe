compfactor <- function(fl, WW=FALSE) {
  purefls <- sapply(fl,function(f) is.null(attr(f,'x')))
  fl <- fl[purefls]
  if(length(fl) <= 1) return(factor(rep(1,length(fl[[1]]))))
  if(WW && length(fl) > 2) {
    cf <- factor(.Call(C_wwcomp,fl))
  } else {
    cf <- factor(.Call(C_conncomp,fl[1:2]))
  }
  cf
}
