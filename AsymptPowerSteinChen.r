#Computing asymptotic power curves using Theorem 3.1.


#Functions for ZIP distribution

#Default parametrization (lambda,omega):

dzip <- function(x, lambda, omega=0) omega*(x==0) + (1-omega)*dpois(x, lambda)

pzip <- function(x, lambda, omega=0) omega*(x>=0) + (1-omega)*ppois(x, lambda)

rzip <- function(n, lambda, omega=0) rbinom(n, 1, 1-omega)*rpois(n, lambda)

qzip <- function(ql, lambda, omega=0){
	up <- qpois(ql, lambda) #ZIP quantiles not later than Poi quantiles
	cdf <- pzib(0:up, lambda, omega)
	min((0:up)[cdf>=ql])
}
qzip <- Vectorize(qzip)



set.seed(123)
test <- rpois(10, 2.5)

#Sample Stein-Chen statistic:
sTSC <- function(x, f){
	mean(x*f(x)) / mean(x) / mean(f(x+1))
}
sTSC(test, function(x) exp(-x)) #0.8258065


#Approximate moment computation:
mu.f <- function(klvec, dist, f){
	upper <- length(dist)-1
	x <- 0:upper
	
	k <- klvec[1]
	
	mom <- NA
	if(length(klvec)==1) mom <- sum(x^k * dist)
	if(length(klvec)>1){
		l <- klvec[-1]
		mom <- sum( apply(cbind(x^k, apply(outer(x, l, "+"), 2, f)), 1, prod) * dist )
	}
	mom
}
mu.f(c(0,1), dpois(0:10, 0.994), function(x) exp(-x)) #0.1962575
mu.f(c(2), dpois(0:10, 0.994), function(x) exp(-x)) #1.982035

#True Stein-Chen statistic:
TSC <- function(dist, f){
	mu.f(c(1,0), dist, f) / mu.f(1, dist, f) / mu.f(c(0,1), dist, f)
}
TSC(dpois(0:10, 0.994), function(x) exp(-x)) #1
TSC(dnbinom(0:50, 10, 0.4), function(x) exp(-x)) #0.513345



#Asymptotics of test statistics:

#"Normalized" covariances, i.e., Cov/mu1/mu2:
nsig <- function(dist, f){
	
	#Some expressions occur repeatedly:
	mu <- mu.f(1, dist, f)
	mu10 <- mu.f(c(1,0), dist, f)
	mu01 <- mu.f(c(0,1), dist, f)

	sig11 <- mu.f(2, dist, f)/mu^2 - 1
	sig12 <- mu.f(c(2,0), dist, f)/mu/mu10 - 1
	sig13 <- mu.f(c(1,1), dist, f)/mu/mu01 - 1
	sig22 <- mu.f(c(2,0,0), dist, f)/mu10^2 - 1
	sig23 <- mu.f(c(1,0,1), dist, f)/mu10/mu01 - 1
	sig33 <- mu.f(c(0,1,1), dist, f)/mu01^2 - 1
	
	c(sig11, sig12, sig13, sig22, sig23, sig33)
}

#Asymptotic mean and SD:
TSCmeansd <- function(n, dist, f){
	Tf <- TSC(dist, f)
	nsigs <- nsig(dist, f)
	
	muTf <- Tf + Tf/n * sum(nsigs * c(1,-1,1,0,-1,1))
	sigTf <- sqrt(Tf^2/n * sum(nsigs * c(1,-2,2,1,-2,1)))
	
	c(muTf, sigTf)
}
TSCmeansd(156, dpois(0:10, 0.994), function(x) exp(-x)) #1.00602803 0.08339678
TSCmeansd(156, dnbinom(0:50, 10, 0.4), function(x) exp(-x)) #0.8570471 0.3898816





#Testing H0: Poi with mu
#against H1: another dist

#Two-sided critical values:
crit.poi <- function(mu, n, f, upper=50, level=0.05){
	z2 <- qnorm(1-level/2)
	pmf <- dpois(0:upper, mu)
	
	musig <- TSCmeansd(n, pmf, f)
	musig[1] + c(-1,1) * z2 * musig[2]
}
crit.poi(0.994, 156, function(x) exp(-x)) #0.8425732 1.1694827

#Exceedance probability for given interval "crit" under alternative:
apower <- function(n, dist, f, crit){
	musig <- TSCmeansd(n, dist, f)

	pnorm(crit[1], musig[1], musig[2]) + 1-pnorm(crit[2], musig[1], musig[2])
}
apower(156, dnbinom(0:50, 10, 10/10.994), function(x) exp(-x), crit.poi(0.994, 156, function(x) exp(-x))) #0.1067408





#############
#Power curves
#############



#Power with f(x)=|x-1|^a against NB alternative
#note that NB-dispersion index equals I=1+mu/n,
#i.e., n=mu/(I-1)

upper <- 50

mu <- 2
tabI <- c(1.1, 1.25, 1.5)
tabN <- mu/(tabI-1)
no.I <- length(tabI)
# n <- 250
n <- 100

tab.a <- seq(0.1,3,0.1)
no.a <- length(tab.a)

tab.crit <- array(NA, c(no.a, 2))
for(k in 1:no.a){
	a <- tab.a[k]
	f <- function(x) abs(x-1)^a
	tab.crit[k,] <- crit.poi(mu, n, f)
}
#tab.crit[1:3,]
#n=250:
          # [,1]     [,2]
# [1,] 0.9391202 1.060482
# [2,] 0.9400753 1.059141
# [3,] 0.9397235 1.059115

tab.power <- array(NA, c(no.a, no.I))
for(k in 1:no.a){
	a <- tab.a[k]
	f <- function(x) abs(x-1)^a
	crit <- tab.crit[k,]
	
for(l in 1:no.I){
	N <- tabN[l]
	pmf <- dnbinom(0:upper, N, N/(N+mu))
	tab.power[k,l] <- apower(n, pmf, f, crit)
}}

cols <- grey(c(0,0.3,0.6))

matplot(tab.a, tab.power, type="l", lwd=2, ylim=c(0,1), col=cols, lty=1, xlab="a", ylab="Power", main=paste("H0: Poi; H1: NB; mu=", mu, ", n=", n, sep=""))
legend("topright",legend=tabI, title=expression(I[Poi]*"="), lwd=2, col=cols, lty=1, cex=0.8, bg=grey(1))
abline(h=0.05, lty=2)


#Find optimal a numerically:
for(l in 1:no.I){
	N <- tabN[l]
	pmf <- dnbinom(0:upper, N, N/(N+mu))

	f.optim <- function(a){
		f <- function(x) abs(x-1)^a
		-apower(n, pmf, f, crit.poi(mu, n, f))
	}
	
	a.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(4.9999), control=list(ndeps=c(1e-4)))
	print(c(n,mu,tabI[l], a.opt$par, -a.opt$value))
} #for I
#100   2.0000000   1.1000000   1.7193634   0.1437202
#100   2.0000000   1.2500000   1.0337769   0.4365804
#100   2.0000000   1.5000000   0.8167225   0.8226934

#250   2.0000000   1.1000000   1.1260295   0.2296826
#250   2.0000000   1.2500000   0.9272283   0.7321826
#250   2.0000000   1.5000000   0.7628869   0.9861976




#Power with f(x)=|x-1|^a against ZIP alternative
#note that mu=(1-om)*lam and I=1+om*lam,
#i.e., lam=mu+I-1

upper <- 50

mu <- 2
tabI <- c(1.1, 1.25, 1.5)
tab.lam <- mu+tabI-1
tab.om <- 1-mu/tab.lam
no.I <- length(tabI)
# n <- 250
n <- 100

tab.a <- seq(0.1,3,0.1)
no.a <- length(tab.a)

tab.power <- array(NA, c(no.a, no.I))
for(k in 1:no.a){
	a <- tab.a[k]
	f <- function(x) abs(x-1)^a
	crit <- crit.poi(mu, n, f)
	
for(l in 1:no.I){
	lam <- tab.lam[l]
	om <- tab.om[l]
	pmf <- dzip(0:upper, lam, om)
	tab.power[k,l] <- apower(n, pmf, f, crit)
}}

cols <- grey(c(0,0.3,0.6))

matplot(tab.a, tab.power, type="l", lwd=2, ylim=c(0,1), col=cols, lty=1, xlab="a", ylab="Power", main=paste("H0: Poi; H1: ZIP; mu=", mu, ", n=", n, sep=""))
legend("topright",legend=tabI, title=expression(I[Poi]*"="), lwd=2, col=cols, lty=1, cex=0.8)
abline(h=0.05, lty=2)


#Find optimal a numerically:
for(l in 1:no.I){
	lam <- tab.lam[l]
	om <- tab.om[l]
	pmf <- dzip(0:upper, lam, om)

	f.optim <- function(a){
		f <- function(x) abs(x-1)^a
		-apower(n, pmf, f, crit.poi(mu, n, f))
	}
	
	a.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(4.9999), control=list(ndeps=c(1e-4)))
	print(c(n,mu,tabI[l], a.opt$par, -a.opt$value))
} #for I
#100   2.0000000   1.1000000   0.2410169   0.2113141
#100   2.0000000   1.2500000   0.2312942   0.6976349
#100   2.0000000   1.5000000   0.2172152   0.9782204

#250   2.0000000   1.1000000   0.2376410   0.3927448
#250   2.0000000   1.2500000   0.2281356   0.9529368
#250   2.0000000   1.5000000   0.2133844   0.9999413



#Power with Box-Cox f(x)=(x^b-1)/b against ZIP alternative
#note that mu=(1-om)*lam and I=1+om*lam,
#i.e., lam=mu+I-1

upper <- 50

mu <- 2
tabI <- c(1.1, 1.25, 1.5)
tab.lam <- mu+tabI-1
tab.om <- 1-mu/tab.lam
no.I <- length(tabI)
# n <- 250
n <- 100

tab.b <- seq(0.1,3,0.1)
no.b <- length(tab.b)

tab.power <- array(NA, c(no.b, no.I))
for(k in 1:no.b){
	b <- tab.b[k]
	f <- function(x) (x^b-1)/b
	crit <- crit.poi(mu, n, f)
	
for(l in 1:no.I){
	lam <- tab.lam[l]
	om <- tab.om[l]
	pmf <- dzip(0:upper, lam, om)
	tab.power[k,l] <- apower(n, pmf, f, crit)
}}

cols <- grey(c(0,0.3,0.6))

matplot(tab.b, tab.power, type="l", lwd=2, ylim=c(0,1), col=cols, lty=1, xlab="b", ylab="Power", main=paste("H0: Poi; H1: ZIP; mu=", mu, ", n=", n, sep=""))
legend("topright",legend=tabI, title=expression(I[Poi]*"="), lwd=2, col=cols, lty=1, cex=0.8, bg=grey(1))
abline(h=0.05, lty=2)


#Find optimal b numerically:
for(l in 1:no.I){
	lam <- tab.lam[l]
	om <- tab.om[l]
	pmf <- dzip(0:upper, lam, om)

	f.optim <- function(b){
		f <- function(x) (x^b-1)/b
		-apower(n, pmf, f, crit.poi(mu, n, f))
	}
	
	b.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(4.9999), control=list(ndeps=c(1e-4)))
	print(c(n,mu,tabI[l], b.opt$par, -b.opt$value))
} #for I
#100   2.0000000   1.1000000   0.0001000   0.1811338
#100   2.0000000   1.2500000   0.0001000   0.6202016
#100   2.0000000   1.5000000   0.0001000   0.9597187

#250   2.0000000   1.1000000   0.0001000   0.3311144
#250   2.000000   1.250000   0.000100   0.914075
#250   2.0000000   1.5000000   0.0001000   0.9997307




#Power with Box-Cox f(x)=(x^b-1)/b against NB alternative
#note that NB-dispersion index equals I=1+mu/n,
#i.e., n=mu/(I-1)

upper <- 50

mu <- 2
tabI <- c(1.1, 1.25, 1.5)
tabN <- mu/(tabI-1)
no.I <- length(tabI)
# n <- 250
n <- 100

tab.b <- seq(0.1,3,0.1)
no.b <- length(tab.b)

tab.power <- array(NA, c(no.b, no.I))
for(k in 1:no.b){
	b <- tab.b[k]
	f <- function(x) (x^b-1)/b
	crit <- crit.poi(mu, n, f)
	
for(l in 1:no.I){
	N <- tabN[l]
	pmf <- dnbinom(0:upper, N, N/(N+mu))
	tab.power[k,l] <- apower(n, pmf, f, crit)
}}

cols <- grey(c(0,0.3,0.6))

matplot(tab.b, tab.power, type="l", lwd=2, ylim=c(0,1), col=cols, lty=1, xlab="b", ylab="Power", main=paste("H0: Poi; H1: NB; mu=", mu, ", n=", n, sep=""))
legend("topright",legend=tabI, title=expression(I[Poi]*"="), lwd=2, col=cols, lty=1, cex=0.8, bg=grey(1))
abline(h=0.05, lty=2)


#Find optimal b numerically:
for(l in 1:no.I){
	N <- tabN[l]
	pmf <- dnbinom(0:upper, N, N/(N+mu))

	f.optim <- function(b){
		f <- function(x) (x^b-1)/b
		-apower(n, pmf, f, crit.poi(mu, n, f))
	}
	
	b.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(4.9999), control=list(ndeps=c(1e-4)))
	print(c(n,mu,tabI[l], b.opt$par, -b.opt$value))
} #for I
#100   2.000000   1.100000   4.999900   0.164935
#100   2.0000000   1.2500000   1.1801354   0.4374394
#100   2.0000000   1.5000000   0.6398454   0.8231910

#250   2.000000   1.100000   1.359233   0.230879
#250   2.0000000   1.2500000   0.8841165   0.7319552
#250   2.0000000   1.5000000   0.4914568   0.9868142




#Power with f(x)=u^x against ZIP alternative
#note that mu=(1-om)*lam and I=1+om*lam,
#i.e., lam=mu+I-1

upper <- 50

mu <- 2
tabI <- c(1.1, 1.25, 1.5)
tab.lam <- mu+tabI-1
tab.om <- 1-mu/tab.lam
no.I <- length(tabI)
n <- 100

tab.u <- seq(0.01,0.99,0.01)
no.u <- length(tab.u)

tab.power <- array(NA, c(no.u, no.I))
for(k in 1:no.u){
	u <- tab.u[k]
	f <- function(x) u^x
	crit <- crit.poi(mu, n, f)
	
for(l in 1:no.I){
	lam <- tab.lam[l]
	om <- tab.om[l]
	pmf <- dzip(0:upper, lam, om)
	tab.power[k,l] <- apower(n, pmf, f, crit)
}}

cols <- grey(c(0,0.3,0.6))

matplot(tab.u, tab.power, type="l", lwd=2, ylim=c(0,1), col=cols, lty=1, xlab="u", ylab="Power", main=paste("H0: Poi; H1: ZIP; mu=", mu, ", n=", n, sep=""))
legend("topright",legend=tabI, title=expression(I[Poi]*"="), lwd=2, col=cols, lty=1, cex=0.8, bg=grey(1))
abline(h=0.05, lty=2)


#Find optimal u numerically:
for(l in 1:no.I){
	lam <- tab.lam[l]
	om <- tab.om[l]
	pmf <- dzip(0:upper, lam, om)

	f.optim <- function(u){
		f <- function(x) u^x
		-apower(n, pmf, f, crit.poi(mu, n, f))
	}
	
	u.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(0.9999), control=list(ndeps=c(1e-4)))
	print(c(n,mu,tabI[l], u.opt$par, -u.opt$value))
} #for I
#100   2.0000000   1.1000000   0.5856217   0.1546737
#100   2.0000000   1.2500000   0.4803183   0.5975422
#100   2.0000000   1.5000000   0.3208073   0.9927516




#Power with f(x)=u^x against NB alternative
#note that NB-dispersion index equals I=1+mu/n,
#i.e., n=mu/(I-1)

upper <- 50

mu <- 2
tabI <- c(1.1, 1.25, 1.5)
tabN <- mu/(tabI-1)
no.I <- length(tabI)
n <- 100

tab.u <- seq(0.01,0.99,0.01)
no.u <- length(tab.u)

tab.power <- array(NA, c(no.u, no.I))
for(k in 1:no.u){
	u <- tab.u[k]
	f <- function(x) u^x
	crit <- crit.poi(mu, n, f)
	
for(l in 1:no.I){
	N <- tabN[l]
	pmf <- dnbinom(0:upper, N, N/(N+mu))
	tab.power[k,l] <- apower(n, pmf, f, crit)
}}

cols <- grey(c(0,0.3,0.6))

matplot(tab.u, tab.power, type="l", lwd=2, ylim=c(0,1), col=cols, lty=1, xlab="u", ylab="Power", main=paste("H0: Poi; H1: NB; mu=", mu, ", n=", n, sep=""))
legend("topleft",legend=tabI, title=expression(I[Poi]*"="), lwd=2, col=cols, lty=1, cex=0.8)
abline(h=0.05, lty=2)


#Find optimal u numerically:
for(l in 1:no.I){
	N <- tabN[l]
	pmf <- dnbinom(0:upper, N, N/(N+mu))

	f.optim <- function(u){
		f <- function(x) u^x
		-apower(n, pmf, f, crit.poi(mu, n, f))
	}
	
	u.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(0.9999), control=list(ndeps=c(1e-4)))
	print(c(n,mu,tabI[l], u.opt$par, -u.opt$value))
} #for I
#100   2.0000000   1.1000000   0.9999000   0.1406306
#100   2.0000000   1.2500000   0.9999000   0.4369378
#100   2.0000000   1.5000000   0.9330555   0.8166315








############
#Simulations
############

set.seed(123)
reps <- 1e5

upper <- 50

mu <- 4
tabI <- c(1.1, 1.25, 1.5, 2)
no.I <- length(tabI)
tab.lam <- mu+tabI-1
tab.om <- 1-mu/tab.lam
tabn <- c(50,100,250)
no.n <- length(tabn)



#Finite-sample performance of sizes:
reject <- array(0, c(no.I, no.n, 2))
opt.des <- array(0, c(no.I, no.n, 2))
for(l in 1:no.I){
	Id <- tabI[l]
	lam <- tab.lam[l]
	om <- tab.om[l]
	pmf <- dzip(0:upper, lam, om)

for(k in 1:no.n){
	n <- tabn[k]
	f.optim <- function(u){
		f <- function(x) u^x
		-apower(n, pmf, f, crit.poi(mu, n, f))
	}
	
	u.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(0.9999), control=list(ndeps=c(1e-4)))
	apow <- -u.opt$value
	u.opt <- u.opt$par
	opt.des[l,k,1] <- u.opt
	opt.des[l,k,2] <- apow

	f.opt <- function(x) u.opt^x
	crit.opt <- crit.poi(mu, n, f.opt)
	
	for(r in 1:reps){
		data0 <- rpois(n, mu)
		data1 <- rzip(n, lam, om)
		
		T0 <- sTSC(data0, f.opt)
		T1 <- sTSC(data1, f.opt)
		
		if(T0<crit.opt[1] || T0>crit.opt[2]) reject[l,k,1] <- reject[l,k,1]+1
		if(T1<crit.opt[1] || T1>crit.opt[2]) reject[l,k,2] <- reject[l,k,2]+1
	} #for reps

	#print(c(mu,Id,n, opt.des[l,k,], reject[l,k,]))
}} #for n,I

results <- array(NA, c(no.I*no.n, 7))
i <- 0
for(l in 1:no.I){
for(k in 1:no.n){
	i <- i+1
	results[i,] <- c(mu, tabI[l], tabn[k], opt.des[l,k,], reject[l,k,])
}} #for n,I
results
 # [1,]    2 1.10   50 0.7104067 0.1044773 4812  10498
 # [2,]    2 1.10  100 0.5856217 0.1546737 5024  15585
 # [3,]    2 1.10  250 0.4645754 0.3197718 5091  32598
 # [4,]    2 1.25   50 0.6021724 0.3386218 4964  34293
 # [5,]    2 1.25  100 0.4803183 0.5975422 5056  60741
 # [6,]    2 1.25  250 0.3294196 0.9556228 5117  94836
 # [7,]    2 1.50   50 0.4600201 0.8272551 5342  82203
 # [8,]    2 1.50  100 0.3208073 0.9927516 5439  98372
 # [9,]    2 1.50  250 0.2535488 1.0000000 5241 100000
# [10,]    2 2.00   50 0.3020811 0.9999863 5846  99873
# [11,]    2 2.00  100 0.5490055 1.0000000 4902 100000
# [12,]    2 2.00  250 0.9998993 1.0000000 4857 100000

 # [1,]    4 1.10   50 0.5705571 0.1729977 4591  17882
 # [2,]    4 1.10  100 0.4992224 0.2769540 4820  29395
 # [3,]    4 1.10  250 0.4006673 0.5754064 5036  59352
 # [4,]    4 1.25   50 0.5192182 0.5197197 4545  55162
 # [5,]    4 1.25  100 0.3997381 0.8325658 4915  82673
 # [6,]    4 1.25  250 0.2596039 0.9994418 5429  99293
 # [7,]    4 1.50   50 0.3926924 0.9623563 4263  92272
 # [8,]    4 1.50  100 0.2980231 0.9999861 5729  99615
 # [9,]    4 1.50  250 0.5535780 1.0000000 4822 100000
# [10,]    4 2.00   50 0.4527591 1.0000000 4333  99930
# [11,]    4 2.00  100 0.6647123 1.0000000 4720 100000
# [12,]    4 2.00  250 0.9998864 0.9999997 4829 100000





#Now NB with power weighting:

set.seed(123)
reps <- 1e5

upper <- 50

mu <- 4
tabI <- c(1.1, 1.25, 1.5, 2)
no.I <- length(tabI)
tabN <- mu/(tabI-1)
tabn <- c(50,100,250)
no.n <- length(tabn)



#Finite-sample performance of sizes:
reject <- array(0, c(no.I, no.n, 2))
opt.des <- array(0, c(no.I, no.n, 2))
for(l in 1:no.I){
	Id <- tabI[l]
	N <- tabN[l]
	pmf <- dnbinom(0:upper, N, N/(N+mu))

for(k in 1:no.n){
	n <- tabn[k]
	f.optim <- function(a){
		f <- function(x) abs(x-1)^a
		-apower(n, pmf, f, crit.poi(mu, n, f))
	}
	
	a.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(4.9999), control=list(ndeps=c(1e-4)))
	apow <- -a.opt$value
	a.opt <- a.opt$par
	opt.des[l,k,1] <- a.opt
	opt.des[l,k,2] <- apow

	f.opt <- function(x) abs(x-1)^a.opt
	crit.opt <- crit.poi(mu, n, f.opt)
	
	for(r in 1:reps){
		data0 <- rpois(n, mu)
		data1 <- rnbinom(n, N, N/(N+mu))
		
		T0 <- sTSC(data0, f.opt)
		T1 <- sTSC(data1, f.opt)
		
		if(T0<crit.opt[1] || T0>crit.opt[2]) reject[l,k,1] <- reject[l,k,1]+1
		if(T1<crit.opt[1] || T1>crit.opt[2]) reject[l,k,2] <- reject[l,k,2]+1
	} #for reps

	#print(c(mu,Id,n, opt.des[l,k,], reject[l,k,]))
}} #for n,I

results <- array(NA, c(no.I*no.n, 7))
i <- 0
for(l in 1:no.I){
for(k in 1:no.n){
	i <- i+1
	results[i,] <- c(mu, tabI[l], tabn[k], opt.des[l,k,], reject[l,k,])
}} #for n,I
results

 # [1,]    2 1.10   50 4.9999000 0.1617505  456   1247
 # [2,]    2 1.10  100 1.7193634 0.1437202 4413  12314
 # [3,]    2 1.10  250 1.1260295 0.2296826 4790  21810
 # [4,]    2 1.25   50 1.2865789 0.2968208 4600  24672
 # [5,]    2 1.25  100 1.0337769 0.4365804 4911  40180
 # [6,]    2 1.25  250 0.9272283 0.7321826 5015  72612
 # [7,]    2 1.50   50 0.8902972 0.6125953 5000  57960
 # [8,]    2 1.50  100 0.8167225 0.8226934 4942  82414
 # [9,]    2 1.50  250 0.7628869 0.9861976 5177  99176
# [10,]    2 2.00   50 0.6521596 0.8978719 5264  92040
# [11,]    2 2.00  100 0.6187068 0.9864418 5132  99509
# [12,]    2 2.00  250 0.5890064 0.9999768 5052 100000

 # [1,]    4 1.10   50 4.9999000 0.1261465 1260   2737
 # [2,]    4 1.10  100 1.4502447 0.1382852 4788  12978
 # [3,]    4 1.10  250 1.0961717 0.2264286 4858  21722
 # [4,]    4 1.25   50 1.1906600 0.2897198 4614  25529
 # [5,]    4 1.25  100 0.9973697 0.4356361 4847  40905
 # [6,]    4 1.25  250 0.9137634 0.7386603 5049  73142
 # [7,]    4 1.50   50 0.8451864 0.6212850 4814  59480
 # [8,]    4 1.50  100 0.8033865 0.8353214 4901  83777
 # [9,]    4 1.50  250 0.7747196 0.9891988 5004  99324
# [10,]    4 2.00   50 0.6571536 0.9140085 4829  93482
# [11,]    4 2.00  100 0.6500828 0.9905665 4779  99664
# [12,]    4 2.00  250 0.6423011 0.9999910 4860 100000

