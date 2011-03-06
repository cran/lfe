library(lfe)
options(lfe.threads=2)
set.seed(655320)
x <- rnorm(20000,mean=2000)
x2 <- rnorm(length(x))
x3 <- rexp(length(x))
## create individual and firm
id <- factor(sample(15000,length(x),replace=TRUE))
firm <- factor(sample(13000,length(x),replace=TRUE))

## effects
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))

## left hand side
y <- x + 0.25*x2 + 0.5*x3 + id.eff[id] + firm.eff[firm] + rnorm(length(x))

## estimate
summary(est <- felm(y ~ x+x2 + x3 + G(id) + G(firm)))
cat('Components:',nlevels(est$cfactor),'largest:',sum(est$cfactor == '1'),'\n')
## extract the group fixed effects
fe <- getfe(est)
## merge back

ideff <- fe[paste('id',id,sep='.'),'effect']
firmeff <- fe[paste('firm',firm,sep='.'),'effect']

## verify that id and firm coefficients are 1
options(scipen=8)
summary(lm(y ~ x + x2 + x3 + ideff + firmeff-1),digits=8)
