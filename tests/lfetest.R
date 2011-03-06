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

## estimate and print result
est <- felm(y ~ x+x2,fl=list(id=id,firm=firm),data=dset)
print(est)

## extract the group fixed effects, and print them
print(getfe(est))
