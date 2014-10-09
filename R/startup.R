.onLoad <- function(libname,pkgname) {
  
  if(is.null(getOption('lfe.usecg')))
    options(lfe.usecg=FALSE)
  if(is.null(getOption('lfe.eps')))
    options(lfe.eps=1e-8)
  if(is.null(getOption('lfe.pint')))
    options(lfe.pint=300)
  if(is.null(getOption('lfe.accel')))
    options(lfe.accel=1)
  if(is.null(getOption('lfe.bootmem')))
    options(lfe.bootmem=500)
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

# acml messes with the affinity, we should store what it is at startup
# on second thought, amd should fix it. A sysadmin may restrict our
# affinity for whatever useful purpose he may have, and we shouldn't
# mess it back (it's enough with acml)
# if(is.null(getOption('lfe.cpuaffinity')))
# options(lfe.cpuaffinity=parallel::mcaffinity())

}

.onUnload <- function(libpath) {
  options(lfe.usecg=NULL, lfe.eps=NULL,lfe.pint=NULL,lfe.accel=NULL,lfe.bootmem=NULL,lfe.threads=NULL,lfe.cpuaffinity=NULL)
  library.dynam.unload('lfe',libpath)
}

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
