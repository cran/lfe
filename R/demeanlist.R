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

