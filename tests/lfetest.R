library(lfe)
options(lfe.threads=4)
set.seed(655318)
x <- rnorm(500)
x2 <- rnorm(length(x))

## create individual and firm
id <- factor(sample(400,length(x),replace=TRUE))
firm <- factor(sample(300,length(x),replace=TRUE,prob=c(2,rep(1,299))))

## effects
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

## left hand side
y <- x + 0.25*x2 + id.eff[id] + firm.eff[firm] + rnorm(length(x))

## estimate and print result
summary(est <- felm(y ~ x+x2+G(id)+G(firm)))

## extract the group fixed effects
getfe(est)
