demeanlist <- function(mtx,fl,icpt=0,eps=getOption('lfe.eps'),
                       threads=getOption('lfe.threads'),
		       progress=getOption('lfe.pint'),
                       accel=getOption('lfe.accel'),
                       randfact=TRUE,
                       means=FALSE) {

  if(length(fl) == 0) {
    if(means) {
      foo <- unlist(as.relistable(mtx))
      foo[] <- 0
      return(relist(foo))
    }
    return(mtx)
  }
  if(is.null(threads)) threads <- 1
  islist <- is.list(mtx)
  if(!islist) mtx <- list(mtx)
# randomize factor order, this may improve convergence
  if(randfact && length(fl) > 2) fl <- fl[order(runif(length(fl)))]

  res <- .Call(C_demeanlist,
     mtx,
     as.list(fl),
     as.integer(icpt),               
     as.double(eps),
     as.integer(threads),
     as.integer(progress),
     as.integer(accel),
     as.logical(means))

  if(!islist) {
    res <- res[[1]]
    names(res) <- names(mtx)
  }
  res
}

# A version of demeanlist with variable argument list. Its only
# purpuse is to avoid a copy of large vectors/matrices when using demeanlist
# on multiple arguments (which is done to enable parallelization over all vectors).
# I.e. demeanlist(list(y,modelmatrix),...) makes a copy of both y and modelmatrix,
# wheras edemeanlist(y,modelmatrix,...) does not.
# Not exported.
edemeanlist <- function(...,fl,icpt=0,eps=getOption('lfe.eps'),
                       threads=getOption('lfe.threads'),
		       progress=getOption('lfe.pint'),
                       accel=getOption('lfe.accel'),
                       randfact=TRUE) {

  if(length(fl) == 0) return(list(...))
  if(is.null(threads)) threads <- 1
# randomize factor order, this may improve convergence
  if(randfact && length(fl) > 2) fl <- fl[order(runif(length(fl)))]

  .External(C_edemeanlist,
     ...,
     fl=as.list(fl),
     icpt=as.integer(icpt),               
     eps=as.double(eps),
     threads=as.integer(threads),
     progress=as.integer(progress),
     accel=as.integer(accel),
     means=FALSE)
}
