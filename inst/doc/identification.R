### R code from vignette source 'identification.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: id (eval = FALSE)
###################################################
## lm(y ~ x1 + x2 + ... + f1 + f2 + ... + fn)


###################################################
### code chunk number 2: id (eval = FALSE)
###################################################
## felm(y ~ x1 + x2 + ... + G(f1) + G(f2) + ... + G(fn))


###################################################
### code chunk number 3: identification.Rnw:83-91
###################################################
library(lfe)
set.seed(42)
x1 <- rnorm(20)
f1 <- sample(8,length(x1),replace=TRUE)/10
f2 <- sample(8,length(x1),replace=TRUE)/10
e1 <- sin(f1) + 0.02*f2^2 + rnorm(length(x1))
y <-  2.5*x1 + (e1-mean(e1))
summary(est <- felm(y ~ x1 + G(f1) + G(f2)))


###################################################
### code chunk number 4: identification.Rnw:95-98
###################################################
ef <- efactory(est,'ref')
is.estimable(ef,est$fe)
getfe(est)


###################################################
### code chunk number 5: identification.Rnw:104-105
###################################################
data.frame(f1,f2,comp=est$cfactor)


###################################################
### code chunk number 6: id (eval = FALSE)
###################################################
## est <- felm(logwage ~ x1 + x2 + G(id) + G(firm) + G(edu))
## getfe(est)


###################################################
### code chunk number 7: id (eval = FALSE)
###################################################
## logwage ~ x1 + x2 + G(id) + G(firm) + edu


###################################################
### code chunk number 8: id (eval = FALSE)
###################################################
## logwage ~ x1 + x2 + G(firm) + G(edu) + G(id)


###################################################
### code chunk number 9: id (eval = FALSE)
###################################################
## y ~ x1 + x2 + f1 + f2 + f3


###################################################
### code chunk number 10: id (eval = FALSE)
###################################################
## est <- felm(y ~ x1 + x2 + G(f1)+G(f2)+G(f3))


###################################################
### code chunk number 11: ex
###################################################
library(lfe)
x1 <- rnorm(100)
f1 <- sample(7,100,replace=TRUE)
f2 <- sample(8,100,replace=TRUE)/8
f3 <- sample(10,100,replace=TRUE)/10
e1 <- sin(f1) + 0.02*f2^2  + 0.17*f3^3 + rnorm(100)
y <-  2.5*x1 + (e1-mean(e1))
summary(est <- felm(y ~ x1 + G(f1) + G(f2) + G(f3)))


###################################################
### code chunk number 12: ex
###################################################
ef <- function(gamma,addnames) {
  ref2 <- gamma[[8]]
  ref3 <- gamma[[16]]
  gamma[1:7] <- gamma[1:7]+ref2+ref3
  gamma[8:15] <- gamma[8:15]-ref2
  gamma[16:25] <- gamma[16:25]-ref3
  if(addnames) {
    names(gamma) <- c(paste('f1',1:7,sep='.'),
                          paste('f2',1:8,sep='.'),
                          paste('f3',1:10,sep='.'))
  }
  gamma
}
is.estimable(ef,fe=est$fe)
getfe(est,ef=ef)


###################################################
### code chunk number 13: ex
###################################################
getfe(est)


###################################################
### code chunk number 14: ex
###################################################
efactory(est,'ref')


###################################################
### code chunk number 15: identification.Rnw:276-284
###################################################
set.seed(42)
x1 <- rnorm(100)
f1 <- sample(34,100,replace=TRUE)
f2 <- sample(34,100,replace=TRUE)/8
f3 <- sample(34,100,replace=TRUE)/10
e1 <- sin(f1) + 0.02*f2^2  + 0.17*f3^3 + rnorm(100)
y <-  2.5*x1 + (e1-mean(e1))
summary(est <- felm(y ~ x1 + G(f1) + G(f2) + G(f3)))


###################################################
### code chunk number 16: identification.Rnw:290-292
###################################################
ef <- efactory(est,'ref')
is.estimable(ef,est$fe)


###################################################
### code chunk number 17: identification.Rnw:300-308
###################################################
f1 <- factor(f1); f2 <- factor(f2); f3 <- factor(f3)
D <- t(do.call('rBind',
        lapply(list(f1,f2,f3),as,Class='sparseMatrix')))
dim(D)
as.integer(rankMatrix(D))
# alternatively we can use an internal function
# in lfe for finding the rank deficiency directly
lfe:::rankDefic(list(f1,f2,f3))


###################################################
### code chunk number 18: identification.Rnw:320-321
###################################################
summary(est <- felm(y ~ x1 + G(f1) + G(f2) + G(f3), exactDOF=TRUE))


###################################################
### code chunk number 19: identification.Rnw:329-331
###################################################
summary(est <- felm(y ~ x1 + G(f2) + G(f3) + f1))
getfe(est)


###################################################
### code chunk number 20: identification.Rnw:339-340
###################################################
data.frame(f1,f2,f3)[1,]


###################################################
### code chunk number 21: identification.Rnw:345-346
###################################################
summary(est <- lm(y ~ x1 + f1 + f2 + f3))


