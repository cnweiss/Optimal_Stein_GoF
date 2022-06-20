#Stein statistic for NB-null.

#Computing asymptotic power curves using Appendix A.5.


#Functions for ZIP distribution with
#default parametrization (lambda,omega):

#Probability mass function in x:
dzip <- function(x, lambda, omega=0) omega*(x==0) + (1-omega)*dpois(x, lambda)

#Cumulative distribution function in x:
pzip <- function(x, lambda, omega=0) omega*(x>=0) + (1-omega)*ppois(x, lambda)

#n i.i.d. random numbers:
rzip <- function(n, lambda, omega=0) rbinom(n, 1, 1-omega)*rpois(n, lambda)

#Quantile function in quantile level ql:
qzip <- function(ql, lambda, omega=0){
	up <- qpois(ql, lambda) #ZIP quantiles not later than Poi quantiles
	cdf <- pzib(0:up, lambda, omega)
	min((0:up)[cdf>=ql])
}
qzip <- Vectorize(qzip)
#These functions are used for the alternative hypothesis.




#Now, Stein statistics and asymptotics from manuscript.

#Sample NB-Stein statistic (A.4), data x, NB-parameter nu, function f:
sTNBS <- function(x, nu, f){
	(nu/mean(x)+1) * mean(x*f(x)) / mean((nu+x)*f(x+1))

}

#Approximate moment computation (10) of order "klvec" and function f, 
#where upper truncation limit M implicitly provided by length of pmf vector "dist":
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
# mu.f(c(0,1), dnbinom(0:10, 1, 1/(1+2.5)), function(x) exp(-x)) #0.1425722
# mu.f(c(2), dnbinom(0:10, 1, 1/(1+2.5)), function(x) exp(-x)) #10.28344

#Population NB-Stein statistic (A.5),
#for pmf vector "dist", NB-parameter nu and function f:
TNBS <- function(dist, nu, f){
	(nu/mu.f(1, dist, f)+1) * mu.f(c(1,0), dist, f) / (nu*mu.f(c(0,1), dist, f) + mu.f(c(1,1), dist, f))
}
# TNBS(dpois(0:10, 2.5), 10, function(x) exp(-x)) #1.144974
# TNBS(dnbinom(0:50, 10, 0.8), 10, function(x) exp(-x)) #1



#Asymptotics of test statistics according to Appendix A.5:

#"Normalized" covariances, i.e., Cov/mu1/mu2,
#for pmf vector "dist", NB-parameter nu and function f:
nsig <- function(dist, nu, f){
	
	#Some expressions occur repeatedly:
	mu <- mu.f(1, dist, f)
	mu10 <- mu.f(c(1,0), dist, f)
	numu01mu11 <- nu*mu.f(c(0,1), dist, f) + mu.f(c(1,1), dist, f)

	sig11 <- mu.f(2, dist, f)/mu^2 - 1
	sig12 <- mu.f(c(2,0), dist, f)/mu/mu10 - 1
	sig22 <- mu.f(c(2,0,0), dist, f)/mu10^2 - 1
	
	#The following differ from the Poisson case:
	sig13 <- (nu*mu.f(c(1,1), dist, f) + mu.f(c(2,1), dist, f)) / mu / numu01mu11 - 1
	sig23 <- (nu*mu.f(c(1,0,1), dist, f) + mu.f(c(2,0,1), dist, f)) / mu10 / numu01mu11 - 1
	sig33 <- (nu^2*mu.f(c(0,1,1), dist, f) + 2*nu*mu.f(c(1,1,1), dist, f) + mu.f(c(2,1,1), dist, f)) / numu01mu11^2 - 1
	
	c(sig11, sig12, sig13, sig22, sig23, sig33)
}

#Asymptotic mean and SD with sample size n,
#for pmf vector "dist", NB-parameter nu and function f:
TNBSmeansd <- function(n, dist, nu, f){
	mu <- mu.f(1, dist, f)

	Tf <- TNBS(dist, nu, f)
	nsigs <- nsig(dist, nu, f)
	
	muTf <- Tf + Tf/n * sum(nsigs * c(nu/(nu+mu), -nu/(nu+mu), nu/(nu+mu), 0,-1,1))
	sigTf <- sqrt(Tf^2/n * sum(nsigs * c(nu^2/(nu+mu)^2, -2*nu/(nu+mu), 2*nu/(nu+mu), 1,-2,1)))
	
	c(muTf, sigTf)
}
# TNBSmeansd(156, dpois(0:10, 2.5), 10, function(x) exp(-x)) #1.1554766 0.1119096
# TNBSmeansd(156, dnbinom(0:50, 10, 0.8), 10, function(x) exp(-x)) #1.0097497 0.1036744





#Testing H0: NB with mu, nu
#against H1: another dist

#Two-sided critical values for given level,
#with mean mu, NB-parameter nu, sample size n, function f, and upper truncation:
crit.nb <- function(mu,nu, n, f, upper=50, level=0.05){
	z2 <- qnorm(1-level/2)
	pmf <- dnbinom(0:upper, nu, nu/(nu+mu))
	
	musig <- TNBSmeansd(n, pmf, nu, f)
	musig[1] + c(-1,1) * z2 * musig[2]
}
#crit.nb(2.5,10, 156, function(x) exp(-x)) #0.8065517 1.2129478

#Exceedance probability for given interval "crit" under alternative (specified by pmf vector "dist"),
#with sample size n, NB-parameter nu and function f,
#where crit might be computed using function "crit.nb":
apower <- function(n, dist, nu, f, crit){
	musig <- TNBSmeansd(n, dist, nu, f)

	pnorm(crit[1], musig[1], musig[2]) + 1-pnorm(crit[2], musig[1], musig[2])
}
# apower(156, dnbinom(0:50, 10, 0.8), 10, function(x) exp(-x), crit.nb(2.5,10, 156, function(x) exp(-x))) #0.05
# apower(156, dpois(0:10, 2.5), 10, function(x) exp(-x), crit.nb(2.5,10, 156, function(x) exp(-x))) #0.304694





##############################
#Power curves of Appendix A.5:
##############################



#Power curves of Figure A1.

#Power with f(x)=|x-1|^a against Poi alternative
#note that NB-dispersion index equals I=1+mu/nu,
#i.e., nu=mu/(I-1)

upper <- 50

mu <- 2
tabI <- c(1.25, 1.5, 2)
tabnu <- mu/(tabI-1)
no.I <- length(tabI)
n <- 100

tab.a <- seq(0.1,3,0.1)
no.a <- length(tab.a)

tab.crit <- array(NA, c(no.a, no.I, 2))
for(l in 1:no.I){
	nu <- tabnu[l]
for(k in 1:no.a){
	a <- tab.a[k]
	f <- function(x) abs(x-1)^a
	tab.crit[k,l,] <- crit.nb(mu,nu, n, f)
}}

tab.power <- array(NA, c(no.a, no.I))
pmf <- dpois(0:upper, mu)
for(l in 1:no.I){
	nu <- tabnu[l]
	
for(k in 1:no.a){
	a <- tab.a[k]
	f <- function(x) abs(x-1)^a
	crit <- tab.crit[k,l,]
	
	tab.power[k,l] <- apower(n, pmf, nu, f, crit)
}}

cols <- grey(c(0,0.3,0.6))

#Figure A1(a):
matplot(tab.a, tab.power, type="l", lwd=2, ylim=c(0,1), col=cols, lty=1, xlab="a", ylab="Power", main=paste("H0: NB; H1: Poi; mu=", mu, ", n=", n, sep=""))
legend("topright",legend=tabI, title=expression(I[Poi]*"="), lwd=2, col=cols, lty=1, cex=0.8, bg=grey(1))
abline(h=0.05, lty=2)


#Find optimal a numerically:
pmf <- dpois(0:upper, mu)
for(l in 1:no.I){
	nu <- tabnu[l]

	f.optim <- function(a){
		f <- function(x) abs(x-1)^a
		-apower(n, pmf, nu, f, crit.nb(mu,nu, n, f))
	}
	
	a.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(4.9999), control=list(ndeps=c(1e-4)))
	print(c(n,mu,tabI[l], a.opt$par, -a.opt$value))
} #for I
#100   2.0000000   1.2500000   0.7373438   0.2284970
#100   2.0000000   1.5000000   0.7374513   0.6752440
#100   2.000000   2.000000   0.734309   0.994752






#Power with f(x)=u^x against ZIP alternative
#note that mu=(1-om)*lam and I=1+om*lam,
#i.e., lam=mu+I-1

upper <- 50

mu <- 2
tabI <- c(1.25, 1.5, 2)
tabnu <- mu/(tabI-1)
tab.lam <- mu+tabI-1
tab.om <- 1-mu/tab.lam
no.I <- length(tabI)
n <- 100

tab.u <- seq(0.01,0.99,0.01)
no.u <- length(tab.u)

tab.power <- array(NA, c(no.u, no.I))
for(l in 1:no.I){
	nu <- tabnu[l]
	lam <- tab.lam[l]
	om <- tab.om[l]
	pmf <- dzip(0:upper, lam, om)
	
for(k in 1:no.u){
	u <- tab.u[k]
	f <- function(x) u^x
	crit <- crit.nb(mu,nu, n, f)
	
	tab.power[k,l] <- apower(n, pmf, nu, f, crit)
}}

cols <- grey(c(0,0.3,0.6))

#Figure 1(b):
matplot(tab.u, tab.power, type="l", lwd=2, ylim=c(0,1), col=cols, lty=1, xlab="u", ylab="Power", main=paste("H0: NB; H1: ZIP; mu=", mu, ", n=", n, sep=""))
legend("topright",legend=tabI, title=expression(I[Poi]*"="), lwd=2, col=cols, lty=1, cex=0.8, bg=grey(1))
abline(h=0.05, lty=2)


#Find optimal u numerically:
for(l in 1:no.I){
	nu <- tabnu[l]
	lam <- tab.lam[l]
	om <- tab.om[l]
	pmf <- dzip(0:upper, lam, om)

	f.optim <- function(u){
		f <- function(x) u^x
		-apower(n, pmf, nu, f, crit.nb(mu,nu, n, f))
	}
	
	u.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(0.9999), control=list(ndeps=c(1e-4)))
	print(c(n,mu,tabI[l], u.opt$par, -u.opt$value))
} #for I
#100   2.0000000   1.2500000   0.2396737   0.1914678
#100   2.0000000   1.5000000   0.1700921   0.5842271
#100   2.0000000   2.0000000   0.1111926   0.9938710


