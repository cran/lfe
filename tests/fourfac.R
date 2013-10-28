library(lfe)
options(lfe.threads=1)
set.seed(655320)
x <- rnorm(50000,mean=200)
x2 <- rnorm(length(x))
x3 <- rexp(length(x))
## create individual and firm
id <- factor(sample(15000,length(x),replace=TRUE))
firm <- factor(sample(13000,length(x),replace=TRUE))
shoe <- factor(sample(1000,length(x),replace=TRUE))
shirt <- factor(sample(1000,length(x),replace=TRUE))
## effects
id.eff <- rnorm(nlevels(id))
firm.eff <- rnorm(nlevels(firm))
shoe.eff <- rnorm(nlevels(shoe))
shirt.eff <- rnorm(nlevels(shirt))
## left hand side
y <- x + 0.25*x2 + 0.5*x3 + id.eff[id] + firm.eff[firm] + shoe.eff[shoe] + shirt.eff[shirt] + rnorm(length(x))

## estimate
summary(est <- felm(y ~ x+x2 + x3 + G(id) + G(firm) + G(shoe) + G(shirt)))
cat('Components:',nlevels(est$cfactor),'largest:',sum(est$cfactor == '1'),'\n')
## extract the group fixed effects
  ## verify that id and firm coefficients are 1
options(scipen=8)

for(ef in c('ln','ref','zm','zm2')) {
  fe <- getfe(est,ef=ef)
  ## merge back

  ideff <- fe[paste('id',id,sep='.'),'effect']
  firmeff <- fe[paste('firm',firm,sep='.'),'effect']
  shoeeff <- fe[paste('shoe',shoe,sep='.'),'effect']
  shirteff <- fe[paste('shirt',shirt,sep='.'),'effect']
  if(ef %in% c('zm','zm2')) {
    icpt <- fe[paste('icpt',1:nlevels(est$cfactor),sep='.'),'effect'][est$cfactor]
    lmres <- lm(y ~ x + x2 + x3 + ideff + firmeff + shoeeff +shirteff + icpt-1)
    acc <- coef(lmres)[c('ideff','firmeff','shoeeff','shirteff','icpt')]
  } else {
    lmres <- lm(y ~ x + x2 + x3 + ideff + firmeff + shoeeff +shirteff-1)
    acc <- coef(lmres)[c('ideff','firmeff','shoeeff','shirteff')]
  }
  print(summary(lmres,digits=8))	
  cat('accuracy:',sprintf('%.8e',acc),'\n')
}

