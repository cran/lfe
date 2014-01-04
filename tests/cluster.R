library(lfe)
# From http://diffuseprior.wordpress.com/2012/06/15/standard-robust-and-clustered-standard-errors-computed-in-r/
set.seed(42)
options(lfe.threads=1, digits=3)
ols <- function(form, data, robust=FALSE, cluster=NULL,digits=3){
    r1 <- lm(form, data)
      if(length(cluster)!=0){
            data <- na.omit(data[,c(colnames(r1$model),cluster)])
                r1 <- lm(form, data)
          }
      X <- model.matrix(r1)
      n <- dim(X)[1]
      k <- dim(X)[2]
      if(robust==FALSE & length(cluster)==0){
            se <- sqrt(diag(solve(crossprod(X)) * as.numeric(crossprod(resid(r1))/(n-k))))
                res <- cbind(coef(r1),se)
          }
      if(robust==TRUE){
            u <- matrix(resid(r1))
                meat1 <- t(X) %*% diag(diag(crossprod(t(u)))) %*% X
                dfc <- n/(n-k)
                se <- sqrt(dfc*diag(solve(crossprod(X)) %*% meat1 %*% solve(crossprod(X))))
                res <- cbind(coef(r1),se)
          }
      if(length(cluster)!=0){
            clus <- cbind(X,data[,cluster],resid(r1))
                colnames(clus)[(dim(clus)[2]-1):dim(clus)[2]] <- c(cluster,"resid")
                m <- dim(table(clus[,cluster]))
                dfc <- (m/(m-1))*((n-1)/(n-k))
                uclust  <- apply(resid(r1)*X,2, function(x) tapply(x, clus[,cluster], sum))
                se <- sqrt(diag(solve(crossprod(X)) %*% (t(uclust) %*% uclust) %*% solve(crossprod(X)))*dfc)
                res <- cbind(coef(r1),se)
          }
      res <- cbind(res,res[,1]/res[,2],(1-pnorm(abs(res[,1]/res[,2])))*2)
      res1 <- matrix(as.numeric(sprintf(paste("%.",paste(digits,"f",sep=""),sep=""),res)),nrow=dim(res)[1])
      rownames(res1) <- rownames(res)
      colnames(res1) <- c("Estimate","Std. Error","t value","Pr(>|t|)")
      return(res1)
  }



x <- rnorm(100)
f1 <- sample(10,length(x), repl=T)
f2 <- sample(10,length(x), repl=T)

y <- x +  log(f1) + cos(f2) + rnorm(length(x), sd=0.5)
dat1 <- data.frame(y, x, f1=factor(f1), f2=factor(f2),cluster=factor(1:length(x)))
#print(summary(lm(y ~ x + f1 + f2, dat1)))
dat2 <- rbind(dat1,dat1)
#summary(lm(y ~ x + f1 + f2, dat2))
#clu <- factor(rep(1:length(x),2))
clu <- dat2[,'cluster']
summary(est <- felm(y ~x + G(f1) + G(f2), dat1))
summary(est <- felm(y ~x + G(f1) + G(f2), dat2, clustervar='cluster'))
summary(est <- felm(y ~x + G(f1) + G(f2), dat2), robust=TRUE)

# then the fixed effects.
ef <- function(gamma, addnames) {
  ref1 <- gamma[[1]]
  ref2 <- gamma[[11]]
  icpt <- ref1+ref2
  res <- c(icpt,gamma[2:10]-ref1,gamma[12:20]-ref2)
  if(addnames) {
    names(res) <- c('icpt',paste('f1',2:10),paste('f2',2:10))
  }
  res
}

ols(y ~x + f1 + f2, dat2, robust=TRUE)
getfe(est,ef=ef,se=T,bN=2000, robust=TRUE)
ols(y ~x + f1 + f2, dat2, cluster="cluster")
getfe(est,ef=ef,se=T,bN=2000, cluster=clu)
