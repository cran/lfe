library(lfe)
options(lfe.threads=4)
set.seed(6553)
# single fixed effect, special case which we risk destroying when we optimize, so check it
x <- rnorm(20000,mean=2000)
x2 <- rnorm(length(x))
x3 <- rexp(length(x))
## create individual and firm
id <- factor(sample(15000,length(x),replace=TRUE))

## effects
id.eff <- rnorm(nlevels(id))

## left hand side
y <- x + 0.25*x2 + 0.5*x3 + id.eff[id] + rnorm(length(x))

## estimate
summary(est <- felm(y ~ x+x2 + x3 + G(id)))

## extract the group fixed effects
fe <- getfe(est)
## merge back

ideff <- fe[paste('id',id,sep='.'),'effect']

## verify that id and firm coefficients are 1
options(scipen=8)
summary(lm(y ~ x + x2 + x3 + ideff -1),digits=8)
