library(lfe)
options(lfe.threads=4)
set.seed(42)
x <- rnorm(1000)
x2 <- rnorm(length(x))

id <- factor(sample(10,length(x),replace=TRUE))
firm <- factor(sample(3,length(x),replace=TRUE,prob=c(2,1,1)))

id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

## left hand side
u <- rnorm(length(x))
x3 <- rnorm(length(x))
Q <- 0.3*x3 + x + 0.2*x2 + id.eff[id] + 0.7*u + rnorm(length(x),sd=0.3)

y <- x + 0.5*x2 + id.eff[id] + firm.eff[firm] + Q + u

## estimate and print result
est <- felm(y ~ x+x2+G(id)+G(firm)+Q,iv=Q~x3)
summary(est,robust=TRUE)

