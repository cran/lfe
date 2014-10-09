# compute the trace of a matrix.
# If we have a list of factors defining a projection, or a function for
# multiplying the matrix with a matrix, we use an iterative method
# from quantum physics, based on the general formula for the expectation
# of a quadratic form E(x' M x) = tr(MV) + E(x)' M E(x) where V=Cov(x). It reduces
# to 
# E(x' M x) = tr(M)
# when x has zero expectation and identity covariance matrix
# the most efficient (lowest variation) is to draw x uniformly from {-1,1}
# "Random phase vector for calculating the trace of a large matrix",
# T. Iitaka and T. Ebisuzaki, Physical Review E 69 (2004).

mctrace <- function(mat, N, tol=1e-3, maxsamples=Inf, ortho=FALSE,
                    trname='',threads=getOption('lfe.threads')) {
  if(is.matrix(mat) || inherits(mat,'Matrix')) {
    return(structure(sum(diag(mat)), sd=0, iterations=0))
  } else if(is.list(mat) && all(sapply(mat, is.factor))) {
    N <- length(mat[[1]])
    if(ortho)
        fun <- function(v) N - colSums(demeanlist(v, mat, threads=threads)*v)
    else
        fun <- function(v) colSums(demeanlist(v, mat, threads=threads)*v)
  } else if(!is.function(mat)) {
    stop('mat must be function, factor list or matrix')
  } else {
    if(missing(N)) stop('N (vector length) must be specified with mat a function')
    if(isTRUE(attr(mat,'IP'))) {
      # inner product is done by the function itself.
      fun <- mat
    } else {
      fun <- function(v) colSums(mat(v)*v)
    }
  }

  if(!is.function(tol)) eps <- function(x) tol else eps <- tol
  
  if(maxsamples < threads) threads <- maxsamples
  maxB <- getOption("lfe.bootmem") * 1e+06
  maxvpt <- maxB %/% (2*8*N*threads)
  if(maxvpt*threads > 4096) maxvpt <- 4096 %/% threads
  
  # ensure at least 4 vectors in first iteration
  if(threads >= 4)
      vpt <- 1
  else
      vpt <- 8 %/% (threads+1)

  blk <- vpt*threads
  if(blk > maxsamples) blk <- maxsamples
  i <- 0
  tr <- 0
  sqsum <- 0
  NN <- 0
  start <- Sys.time()
  last <- 0
  # get a clue about the tolerance.
  cureps <- eps(as.numeric(fun(0, trtol=0)))/2

  while(NN < maxsamples && (NN < 4 || (cureps > 0 && relsd > cureps) || (cureps < 0 && sd > -cureps))) {
    i <- i+1
    now <- Sys.time()
    if(NN > 0) {
      remaining <- as.integer((Ntarget-NN)/(NN/as.numeric(now-start)))
      if(remaining > 60 && now - last > 60 || TRUE) {
        message('  *** trace ',trname,' sample ',NN,' of ',Ntarget,', expected finish at ',
                now + remaining) 
        last <- now
      }
    }
#    if(NN > 0) message('tr=',tr/NN,' sd=',sd,' cureps=',cureps,
#                       ' N=',NN,' blk=',blk,' targ=', Ntarget)
    ests <- fun(matrix(sample(c(-1,1), N*blk, replace=TRUE), N), trtol=abs(cureps))
#    ests <- fun(matrix(rnorm(N*blk), N))
    NN <- NN + blk
    tr <- tr + sum(ests)
    sqsum <- sqsum + sum(ests^2)
    rm(ests)#; gc()
    # compute sd for the mean tr/NN. It's sqrt(E(x^2) - E(X)^2)/sqrt(NN)
    sd <- sqrt(sqsum/NN - (tr/NN)^2)/sqrt(NN)
    if(NN > 1) sd <- sd*sqrt(NN/(NN-1))  # small sample correction?
    relsd <- sd/abs(tr/NN)
    cureps <- eps(tr/NN)
    # try to figure out how many iterations are needed to obtain the
    # desired tolerance.
    sdtarget <- if(cureps < 0) -cureps else cureps*abs(tr/NN)
    Ntarget <- as.integer((sd/sdtarget)^2*NN)
    vpt <- 1 + (Ntarget-NN) %/% threads
    if(vpt > maxvpt) vpt <- maxvpt
    if(vpt < 1) vpt <- 1
    blk <- vpt*threads
  }
#  if(last > start) cat('\n')
  structure(tr/NN, sd=sd, iterations=NN)
}
