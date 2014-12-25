library(lfe)
options(lfe.threads=2,digits=5,warn=1,lfe.eps=1e-5)
set.seed(42)
x <- rnorm(400)
x2 <- rnorm(length(x))

id <- factor(sample(10,length(x),replace=TRUE))
firm <- factor(sample(3,length(x),replace=TRUE,prob=c(2,1,1)))
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

## left hand side
u <- rnorm(length(x))
x3 <- rnorm(length(x))
x4 <- rnorm(length(x))
Q <- 0.3*x3 + x + 0.2*x2 + id.eff[id] + 0.15*u + rnorm(length(x),sd=0.2)
R <- 0.001*x3 + 0.2*x + 0.5*x2 + 0.7*id.eff[id] - 0.11*u + rnorm(length(x),sd=0.2)
y <- x + 0.5*x2 + id.eff[id] + firm.eff[firm] + Q + R + u

## estimate and print result
est <- felm(y ~ x+x2 | id+firm |(Q|R~x3+x4))
summary(est,robust=TRUE)
clu <- factor(sample(10,length(x), replace=TRUE))
# try it from a function
fun <- function() {
  Y <- y
  S <- Q
  felm(Y ~ x+x2 | id + firm |(S|R ~ x3+x4), cluster=clu)
  fr <- data.frame(y,x,x2,id,firm,Q,R,x3,x4)
  felm(y ~ x+x2 | id + firm |(Q|R ~ x3+x4), data=fr, cluster=clu)
}
est <- fun()
summary(est, robust=TRUE)
for(s1 in est$step1) print(summary(s1))
