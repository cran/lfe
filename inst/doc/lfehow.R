### R code from vignette source 'lfehow.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: rmod (eval = FALSE)
###################################################
## Y ~ X1 + X2 + ... + Xk + D1 + D2 + ... + De


###################################################
### code chunk number 2: felm (eval = FALSE)
###################################################
## felm(Y ~ X1 + X2 + ... + Xk + G(D1) + G(D2) + ... + G(De))


###################################################
### code chunk number 3: ex
###################################################
set.seed(41)
x <- rnorm(500)
x2 <- rnorm(length(x))
x3 <- rnorm(length(x))


###################################################
### code chunk number 4: ex
###################################################
f1 <- factor(sample(7,length(x),replace=TRUE))
f2 <- factor(sample(4,length(x),replace=TRUE))
f3 <- factor(sample(3,length(x),replace=TRUE))
eff1 <- rnorm(nlevels(f1))
eff2 <- rexp(nlevels(f2))
eff3 <- runif(nlevels(f3))


###################################################
### code chunk number 5: ex
###################################################
 y <- x + 0.5*x2 + 0.25*x3 + eff1[f1] + eff2[f2] + eff3[f3] + rnorm(length(x))


###################################################
### code chunk number 6: ex
###################################################
demean <- function(v,fl) {
  Pv <- v; oldv <- v-1
  while(sqrt(sum((Pv-oldv)**2)) >= 1e-7) {
    oldv <- Pv
    for(f in fl) Pv <- Pv - ave(Pv,f)
  }
 Pv
}


###################################################
### code chunk number 7: ex
###################################################
 fl <- list(f1,f2,f3)
 Py <- demean(y,fl)
 Px <- demean(x,fl)
 Px2 <- demean(x2,fl)
 Px3 <- demean(x3,fl)


###################################################
### code chunk number 8: ex
###################################################
summary(lm(Py ~ Px + Px2 + Px3 - 1))


###################################################
### code chunk number 9: ex
###################################################
 library(lfe)


###################################################
### code chunk number 10: ex
###################################################
 summary(est <- felm(y ~ x + x2 + x3 + G(f1)+G(f2)+G(f3)))


###################################################
### code chunk number 11: ex
###################################################
ef <- function(v,addnames) {
  r2 <- v[[8]]
  r3 <- v[[12]]
  v[1:7] <- v[1:7] + r2 + r3
  v[8:11] <- v[8:11] - r2
  v[12:14] <- v[12:14] - r3
  if(addnames) names(v) <- c(paste('f1',1:7,sep='.'),
                             paste('f2',1:4,sep='.'),
                             paste('f3',1:3,sep='.'))
  v
}
# verify that it's estimable
is.estimable(ef,est$fe)
getfe(est,ef=ef,se=TRUE)


###################################################
### code chunk number 12: ex
###################################################
summary(lm(y ~ x + x2 + x3 + f1 + f2 + f3 - 1))


