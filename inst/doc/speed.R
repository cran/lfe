
## ------------------------------------------------------------------------
library(lfe)
set.seed(42)
x <- rnorm(100000)
f1 <- sample(10000,length(x),replace=TRUE)
f2 <- sample(10000,length(x),replace=TRUE)
y <- x + cos(f1) + log(f2+1) + rnorm(length(x), sd=0.5)


## ------------------------------------------------------------------------
system.time(est <- felm(y ~ x + G(f1) + G(f2)))


## ------------------------------------------------------------------------
system.time(alpha <- getfe(est))


## ------------------------------------------------------------------------
f2 <- sample(300,length(x),replace=TRUE)
y <- x + cos(f1) + log(f2+1) + rnorm(length(x), sd=0.5)
system.time(est <- felm(y ~ x + G(f1) + G(f2)))
system.time(alpha <- getfe(est))


## ------------------------------------------------------------------------
f2 <- (f1 + sample(20,length(x),replace=TRUE)) %% 300
y <- x + cos(f1) + log(f2+1) + rnorm(length(x), sd=0.5)
system.time(est <- felm(y ~ x + G(f1) + G(f2)))
system.time(alpha <- getfe(est))


## ------------------------------------------------------------------------
f2 <- (f1 + sample(5,length(x),replace=TRUE)) %% 300
y <- x + cos(f1) + log(f2+1) + rnorm(length(x), sd=0.5)
system.time(est <- felm(y ~ x + G(f1) + G(f2)))
system.time(alpha <- getfe(est))
nlevels(est[['cfactor']]) # number of connected components


## ------------------------------------------------------------------------
system.time(est <- felm(y ~ x + G(f1) + factor(f2)))
system.time(alpha <- getfe(est))


## ------------------------------------------------------------------------
f2 <- (f1 + sample(5,length(x),replace=TRUE)^3) %% 300
y <- x + cos(f1) + log(f2+1) + rnorm(length(x), sd=0.5)
system.time(est <- felm(y ~ x + G(f1) + G(f2)))
system.time(alpha <- getfe(est))
nlevels(est[['cfactor']])


## ----cache=TRUE, tidy=FALSE----------------------------------------------
library(igraph)
mkgraph <- function(f1,f2)
  graph.edgelist(cbind(paste('f1',f1),paste('f2',f2)), directed=FALSE)

f2 <- sample(10000,length(x),replace=TRUE)
diameter(mkgraph(f1,f2))
f2 <- (f1 + sample(5,length(x),replace=TRUE)^3) %% 300
diameter(mkgraph(f1,f2))
f2 <- (f1 + sample(5,length(x),replace=TRUE)) %% 300
diameter(mkgraph(f1,f2))


