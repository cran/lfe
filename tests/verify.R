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

## make a data frame
dset <- data.frame(y,x,x2,id,firm)

## estimate
est <- felm(y ~ x+x2,fl=list(id=id,firm=firm),data=dset)

## extract the group fixed effects
fe <- getfe(est)

## merge back into dataset

dset[,'id.eff'] <- fe[paste('id',dset[,'id'],sep='.'),'effect']
dset[,'firm.eff'] <- fe[paste('firm',dset[,'firm'],sep='.'),'effect']

## verify that id and firm coefficients are 1
summary(lm(y ~ x + x2 + id.eff + firm.eff-1,data=dset))
