require('mvtnorm')
# individuals
Nind <- 15000
# firms
Nfirms <- 2000
# and these people change their shoes now and then
# they have 400 different to choose from
# change Nfirms and Nshoes to 1 to simulate stalinist economy
Nshoes <- 400
# avg observations per individual
obs <- 15
years <- 1990:(1990+obs-1)
Nrows <- Nind*obs



#varcorr matrix of observed covariates
vcv <- matrix(c(1,0.5,0.5,0.71),ncol=2)
X <- rmvnorm(Nrows,c(0,0),vcv)
x <- X[,1]
x2 <- X[,2]

# individual effects
ife <- rchisq(Nind,2)
ife <- runif(Nind)

# normalize variance
ife <- ife/sqrt(var(ife))

# firm fixed effects
ffe <- rexp(Nfirms) 
#ffe <- runif(Nfirms)
ffe <- ffe/sqrt(var(ffe))

# shoe fixed effects
sfe <- runif(Nshoes)

yfe <- seq(0.0,(obs-1)/10,length.out=obs)
#yfe <- rnorm(obs)


# hmm, now generate data
# start out with a random assignment to firms
# firms has a probability of being moved to
fprob <- rchisq(Nfirms,10)
sprob <- rchisq(Nshoes,10)
#fprob <- rep(1,Nfirms)
# if want to model assortative matching, we should
# let the probability vary with the fixed effect of
# the individual and the firm.  Hmm, how to do this in a 
# simple fashion?

year = factor(rep(years,Nind))
yint <- as.integer(levels(year))[year]
tab <- data.frame(x=x,x2=x2,year=yint,id=rep(1:Nind,each=obs),firm=0,shoe=0,y=0)
firmy <- sample(1:Nfirms,Nind,replace=T,prob=fprob)
shoey <- sample(1:Nshoes,Nind,replace=T,prob=sprob)
for(i in years) {
   cat('Doing year',i,'\n')
   tab[year == i,'firm'] <- firmy
   tab[year == i,'shoe'] <- shoey
   # Now, figure out who changes  10% prob
   ch <- runif(Nind) < 0.1
   nch <- sum(ch)
   firmy[ch] <- sample(1:Nfirms,nch,replace=T,fprob)
   ch <- runif(Nind) < 0.1
   nch <- sum(ch)
   shoey[ch] <- sample(1:Nshoes,nch,replace=T,sprob)
}

firm <- as.factor(tab[,'firm'])
id <- as.factor(tab[,'id'])
shoe <- as.factor(tab[,'shoe'])
tab[,'ife'] <- ife[id]
tab[,'ffe'] <- ffe[firm]   
tab[,'sfe'] <- sfe[shoe]
tab[,'yfe'] <- yfe[year]
# introduce some correlation
#tab[,'x'] <- tab[,'x'] + 0.3*ife[id] + 0.2*ffe[firm]
#tab[,'x2'] <- tab[,'x2'] + 0.2*ffe[firm] - 0.3*ife[id]
resid <- rnorm(Nrows,sd=sqrt(0.36))
tab[,'y'] <- tab[,'ife'] + tab[,'ffe'] + tab[,'sfe'] + tab[,'yfe'] + 0.5*x + 0.25*x2 + resid
write.table(tab,'shoedata.csv',row.names=FALSE,col.names=FALSE)
if(!interactive()) quit('n')
