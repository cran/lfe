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
### code chunk number 3: identification.Rnw:89-97
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
### code chunk number 4: identification.Rnw:101-104
###################################################
ef <- efactory(est)
is.estimable(ef,est$fe)
getfe(est)


###################################################
### code chunk number 5: identification.Rnw:110-111
###################################################
data.frame(f1,f2,comp=est$cfactor)


###################################################
### code chunk number 6: identification.Rnw:118-120
###################################################
f1 <- factor(f1); f2 <- factor(f2)
summary(lm(y ~ x1 + f1 + f2))


###################################################
### code chunk number 7: id (eval = FALSE)
###################################################
## est <- felm(logwage ~ x1 + x2 + G(id) + G(firm) + G(edu))
## getfe(est)


###################################################
### code chunk number 8: id (eval = FALSE)
###################################################
## logwage ~ x1 + x2 + G(id) + G(firm) + edu


###################################################
### code chunk number 9: id (eval = FALSE)
###################################################
## logwage ~ x1 + x2 + G(firm) + G(edu) + G(id)


###################################################
### code chunk number 10: id (eval = FALSE)
###################################################
## y ~ x1 + x2 + f1 + f2 + f3


###################################################
### code chunk number 11: id (eval = FALSE)
###################################################
## est <- felm(y ~ x1 + x2 + G(f1)+G(f2)+G(f3))


###################################################
### code chunk number 12: ex
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
### code chunk number 13: ex
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
### code chunk number 14: ex
###################################################
getfe(est)


###################################################
### code chunk number 15: ex
###################################################
efactory(est)


###################################################
### code chunk number 16: ex
###################################################
f1 <- factor(f1); f2 <- factor(f2); f3 <- factor(f3)
ef <- function(gamma,addnames) {
  ref1 <- gamma[[1]]
  ref2 <- gamma[[8]]
  ref3 <- gamma[[16]]
  # put the intercept in the first coordinate
  gamma[[1]] <- ref1+ref2+ref3
  gamma[2:7] <- gamma[2:7]-ref1
  gamma[8:14] <- gamma[9:15]-ref2
  gamma[15:23] <- gamma[17:25]-ref3
  length(gamma) <- 23
  if(addnames) {
    names(gamma) <- c('(Intercept)',paste('f1',levels(f1)[2:7],sep=''),
                          paste('f2',levels(f2)[2:8],sep=''),
                          paste('f3',levels(f3)[2:10],sep=''))
  }
  gamma
}
getfe(est,ef=ef,bN=1000,se=TRUE)
#compare with lm
summary(lm(y ~ x1 + f1 + f2 + f3))


###################################################
### code chunk number 17: identification.Rnw:320-328
###################################################
set.seed(128)
x1 <- rnorm(25)
f1 <- sample(9,length(x1),replace=TRUE)
f2 <- sample(8,length(x1),replace=TRUE)
f3 <- sample(8,length(x1),replace=TRUE)
e1 <- sin(f1) + 0.02*f2^2  + 0.17*f3^3 + rnorm(length(x1))
y <-  2.5*x1 + (e1-mean(e1))
summary(est <- felm(y ~ x1 + G(f1) + G(f2) + G(f3)))


###################################################
### code chunk number 18: identification.Rnw:334-336
###################################################
ef <- efactory(est)
is.estimable(ef,est$fe)


###################################################
### code chunk number 19: identification.Rnw:344-349
###################################################
f1 <- factor(f1); f2 <- factor(f2); f3 <- factor(f3)
D <- t(do.call('rBind',
        lapply(list(f1,f2,f3),as,Class='sparseMatrix')))
dim(D)
as.integer(rankMatrix(D))


###################################################
### code chunk number 20: identification.Rnw:353-354
###################################################
lfe:::rankDefic(list(f1,f2,f3))


###################################################
### code chunk number 21: identification.Rnw:366-367
###################################################
summary(est <- felm(y ~ x1 + G(f1) + G(f2) + G(f3), exactDOF=TRUE))


###################################################
### code chunk number 22: identification.Rnw:375-377
###################################################
summary(est <- felm(y ~ x1 + G(f1) + G(f2) + f3, exactDOF=TRUE))
getfe(est)


###################################################
### code chunk number 23: identification.Rnw:385-389
###################################################
summary(est <- felm(y ~ x1 + G(f1) + G(f3) + G(f2), exactDOF=TRUE))
is.estimable(efactory(est),est$fe)
getfe(est)



###################################################
### code chunk number 24: identification.Rnw:399-400
###################################################
data.frame(f1,f2,f3)[1,]


###################################################
### code chunk number 25: identification.Rnw:405-406
###################################################
summary(est <- lm(y ~ x1 + f1 + f2 + f3))


###################################################
### code chunk number 26: identification.Rnw:427-429
###################################################
fe <- list(f1,f2,f3)
wwcomp <- compfactor(fe, WW=TRUE)


###################################################
### code chunk number 27: identification.Rnw:432-434
###################################################
lfe:::rankDefic(fe)
nlevels(wwcomp)


###################################################
### code chunk number 28: identification.Rnw:440-444
###################################################
nlevels(interaction(compfactor(fe),wwcomp))
# pick the largest component:
wwdata <- data.frame(y, x1, f1, f2, f3)[wwcomp==1, ]
print(wwdata)


###################################################
### code chunk number 29: identification.Rnw:465-475
###################################################
set.seed(135)
x <- rnorm(10000)
f1 <- sample(1000,length(x),replace=TRUE)
f2 <- (f1 + sample(18,length(x), replace=TRUE)) %% 500
f3 <- (f2 + sample(9,length(x),replace=TRUE)) %% 500
y <- x + 1e-4*f1 + sin(f2^2) +
  cos(f3)^3 + 0.5*rnorm(length(x))
dataset <- data.frame(y,x,f1,f2,f3)
summary(est <- felm(y ~ x + G(f1) + G(f2) + G(f3),
             data=dataset, exactDOF=TRUE))


###################################################
### code chunk number 30: identification.Rnw:480-483
###################################################
nlevels(est$cfactor)
is.estimable(efactory(est), est$fe)
nrow(alpha <- getfe(est))


###################################################
### code chunk number 31: identification.Rnw:486-487
###################################################
lfe:::rankDefic(est$fe)


###################################################
### code chunk number 32: identification.Rnw:491-497
###################################################
wwcomp <- compfactor(est$fe,WW=TRUE)
nlevels(wwcomp)
wwset <- dataset[wwcomp == 1, ]
nrow(wwset)
summary(wwest <- felm(y ~ x + G(f1) + G(f2) + G(f3),
             data=wwset, exactDOF=TRUE))


###################################################
### code chunk number 33: identification.Rnw:504-505
###################################################
lfe:::rankDefic(wwest$fe)


###################################################
### code chunk number 34: identification.Rnw:508-509
###################################################
nrow(wwalpha <- getfe(wwest))


###################################################
### code chunk number 35: identification.Rnw:513-515
###################################################
head(alpha)
head(wwalpha)


###################################################
### code chunk number 36: identification.Rnw:521-522
###################################################
table(wwalpha[,'fe'])


