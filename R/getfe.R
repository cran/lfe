# return a data-frame with the group fixed effects, including zeros for references
getfe <- function(obj,references=NULL,se=FALSE,method='kaczmarz',ef='ref',bN=100, robust=FALSE, cluster=obj[['clustervar']], lhs=NULL) {

  if(length(obj$fe) == 0) return(NULL)
  if(!is.null(obj$numctrl) && obj$numctrl > 0)
      stop("Can't retrieve fixed effects when estimating with control variables")
  if(method == 'kaczmarz') {
    if(!is.null(references))
       warning('use estimable function (ef) instead of references in the Kaczmarz method')
    if(is.null(ef)) ef <- 'ln'
    if(!is.character(ef) && !is.function(ef))
      stop('ef must be a function when using the Kaczmarz method')
    return(getfe.kaczmarz(obj,se,ef=ef,bN=bN, robust=robust, cluster=cluster, lhs=lhs))
  }
  if(method != 'cholesky') stop('method must be either kaczmarz or cholesky')
  attr(se,'sefactor') <- obj$sefactor
  attr(obj$fe,'references') <- references
  R <- obj$r.residuals
  # then the remaining.  This is usually sufficient.
  # we could also partition differently, just do the 'comp' adjustment accordingly
  # components
  ddlist <- makeddlist(obj$fe)
  gc()
  orignm <- attr(ddlist,'nm')
  comp <- 1
  res <- data.frame()
  for(dd in ddlist) {
#  res <- foreach(dd=ddlist,.combine=rbind,.init=data.frame()) %dopar% {
    dummies <- attr(dd,'dummies')
    keep <- attr(dd,'keep')
    comp <- attr(dd,'comp')
#    cat(date(),'comp dd',comp,'size',length(keep),'\n')
    Rhs <- as.vector(dummies %*% R[keep])
    names(Rhs) <- colnames(dd)
    alpha <- findfe(dd,Rhs,se)
    alpha[,'comp'] <- comp
    res <- rbind(res,alpha)
#    alpha
  }
  res <- res[orignm,]
  res[,'comp'] <- factor(res[,'comp'])    
# now, add factors telling which fe-group we're in
# the rownames are of the form <fe>.<idx>
  fefact <- strsplit(rownames(res),'.',fixed=TRUE)
  res[,'fe'] <- factor(unlist(lapply(fefact,function(l) l[[1]])))
  res[,'idx'] <- factor(unlist(lapply(fefact,function(l) paste(l[-1],collapse='.'))))
  return(res)
}
