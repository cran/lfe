setoption <- function(opt,default) {
  optnam <- paste('lfe',opt,sep='.')
  if(!is.null(getOption(optnam))) return()
  envnam <- paste('LFE',toupper(opt),sep='_')
  e <- Sys.getenv(envnam)
  if(e != '') {
    val <- try(eval(parse(text=e)))
    if(inherits(val,'try-error')) val <- default
  } else val <- default
  ol <- list(val)
  names(ol) <- optnam
  do.call(options,ol)
}

.onLoad <- function(libname,pkgname) {
  
  setoption('usecg',FALSE)
  setoption('eps',1e-8)
  setoption('pint',300)
  setoption('accel',1)
  setoption('bootmem',500)

  if(is.null(getOption('lfe.threads'))) {
    cr <- as.integer(Sys.getenv('LFE_THREADS'))
    if(is.na(cr)) cr <- as.integer(Sys.getenv('OMP_NUM_THREADS'))
    if(is.na(cr)) cr <- as.integer(Sys.getenv('OMP_THREAD_LIMIT'))
    if(is.na(cr)) cr <- as.integer(Sys.getenv('NUMBER_OF_PROCESSORS'))
    if(is.na(cr)) cr <- numcores()
    if(is.na(cr) || cr < 1) {
      cr <- 1
    }
    options(lfe.threads=cr)
  }
}

.onUnload <- function(libpath) {
  options(lfe.usecg=NULL, lfe.eps=NULL,lfe.pint=NULL,lfe.accel=NULL,
          lfe.bootmem=NULL,lfe.threads=NULL)
  library.dynam.unload('lfe',libpath)
}

if(!exists('anyNA')) anyNA <- function(x) any(is.na(x))

numcores <- function() {
  "This function is snatched from package 'multicore'"
  systems <- list(darwin = "/usr/sbin/sysctl -n hw.ncpu 2>/dev/null",
                  freebsd = "/sbin/sysctl -n hw.ncpu 2>/dev/null",
                  linux = "grep processor /proc/cpuinfo 2>/dev/null|wc -l",
                  irix = c("hinv |grep Processors|sed 's: .*::'", "hinv|grep '^Processor '|wc -l"),
                  solaris = "/usr/sbin/psrinfo -v|grep 'Status of.*processor'|wc -l")
  for (i in seq(systems)) if (length(grep(paste("^",
                                                names(systems)[i], sep = ""), R.version$os)))
    for (cmd in systems[i]) {
      a <- gsub("^ +", "", system(cmd, TRUE)[1])
      if (length(grep("^[1-9]", a)))
        return(as.integer(a))
    }
  NA
}
