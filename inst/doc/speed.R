### R code from vignette source 'speed.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: speed.Rnw:52-58
###################################################
library(lfe)
set.seed(42)
x <- rnorm(100000)
f1 <- sample(10000,length(x),replace=TRUE)
f2 <- sample(10000,length(x),replace=TRUE)
y <- x + cos(f1) + log(f2+1) + rnorm(length(x), sd=0.5)


###################################################
### code chunk number 2: speed.Rnw:61-62
###################################################
system.time(est <- felm(y ~ x + G(f1) + G(f2)))


###################################################
### code chunk number 3: speed.Rnw:65-66
###################################################
system.time(alpha <- getfe(est))


###################################################
### code chunk number 4: speed.Rnw:72-76
###################################################
f2 <- sample(300,length(x),replace=TRUE)
y <- x + cos(f1) + log(f2+1) + rnorm(length(x), sd=0.5)
system.time(est <- felm(y ~ x + G(f1) + G(f2)))
system.time(alpha <- getfe(est))


###################################################
### code chunk number 5: speed.Rnw:90-94
###################################################
f2 <- (f1 + sample(20,length(x),replace=TRUE)) %% 300
y <- x + cos(f1) + log(f2+1) + rnorm(length(x), sd=0.5)
system.time(est <- felm(y ~ x + G(f1) + G(f2)))
system.time(alpha <- getfe(est))


###################################################
### code chunk number 6: speed.Rnw:100-105
###################################################
f2 <- (f1 + sample(5,length(x),replace=TRUE)) %% 300
y <- x + cos(f1) + log(f2+1) + rnorm(length(x), sd=0.5)
system.time(est <- felm(y ~ x + G(f1) + G(f2)))
system.time(alpha <- getfe(est))
nlevels(est[['cfactor']]) # number of connected components


###################################################
### code chunk number 7: speed.Rnw:112-114
###################################################
system.time(est <- felm(y ~ x + G(f1) + factor(f2)))
system.time(alpha <- getfe(est))


###################################################
### code chunk number 8: speed.Rnw:125-130
###################################################
f2 <- (f1 + sample(5,length(x),replace=TRUE)^3) %% 300
y <- x + cos(f1) + log(f2+1) + rnorm(length(x), sd=0.5)
system.time(est <- felm(y ~ x + G(f1) + G(f2)))
system.time(alpha <- getfe(est))
nlevels(est[['cfactor']])


###################################################
### code chunk number 9: speed.Rnw:154-167 (eval = FALSE)
###################################################
## mkgraph <- function(f1,f2)
##   graph.adjacency(tcrossprod(rBind(as(factor(f1),'sparseMatrix'),
##                                    as(factor(f2),'sparseMatrix')))>0,
##                   'undirected',diag=FALSE)
## f2 <- sample(10000,length(x),replace=TRUE)
## fivenum(shortest.paths(mkgraph(f1,f2),
##                        v=sample(20000,10),to=sample(20000,10)))
## f2 <- (f1 + sample(5,length(x),replace=TRUE)^3) %% 300
## fivenum(shortest.paths(mkgraph(f1,f2),
##                        v=sample(10300,10), to=sample(10300,10)))
## f2 <- (f1 + sample(5,length(x),replace=TRUE)) %% 300
## fivenum(shortest.paths(mkgraph(f1,f2),
##                        v=sample(10300,10), to=sample(10300,10)))


