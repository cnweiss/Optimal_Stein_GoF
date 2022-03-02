#Computing asymptotic power curves using Theorem 4.1.


#Functions for ZIB distribution

#Default parametrization (p,omega):

dzib <- function(x, N,p, omega=0) omega*(x==0) + (1-omega)*dbinom(x, N,p)

pzib <- function(x, N,p, omega=0) omega*(x>=0) + (1-omega)*pbinom(x, N,p)

rzib <- function(n, N,p, omega=0) rbinom(n, 1, 1-omega)*rbinom(n, N,p)

qzib <- function(ql, N,p, omega=0){
	up <- qbinom(ql, N,p) #ZIB quantiles not later than Bin quantiles
	cdf <- pzib(0:up, N,p, omega)
	min((0:up)[cdf>=ql])
}
qzib <- Vectorize(qzib)






set.seed(123)
test <- rbinom(10, 15,0.15)

#Sample Bin-Stein statistic:
sTBS <- function(x, N, f){
	(N/mean(x)-1) * mean(x*f(x)) / mean((N-x)*f(x+1))
}
sTBS(test, 15, function(x) exp(-x)) #0.8562782


#Approximate moment computation (bounded range):
mu.f <- function(klvec, dist, f){
	N <- length(dist)-1
	x <- 0:N
	
	k <- klvec[1]
	
	mom <- NA
	if(length(klvec)==1) mom <- sum(x^k * dist)
	if(length(klvec)>1){
		l <- klvec[-1]
		mom <- sum( apply(cbind(x^k, apply(outer(x, l, "+"), 2, f)), 1, prod) * dist )
	}
	mom
}
mu.f(c(0,1), dbinom(0:15, 15,0.15), function(x) exp(-x)) #0.08255503
mu.f(c(2), dbinom(0:15, 15,0.15), function(x) exp(-x)) #6.975

#True Bin-Stein statistic:
TBS <- function(dist, f){
	N <- length(dist)-1

	(N/mu.f(1, dist, f)-1) * mu.f(c(1,0), dist, f) / (N*mu.f(c(0,1), dist, f) - mu.f(c(1,1), dist, f))
}
TBS(dbinom(0:15, 15,0.15), function(x) exp(-x)) #1
TBS(dzib(0:15, 15,0.15, 0.4), function(x) exp(-x)) #0.4285461



#Asymptotics of test statistics:

#"Normalized" covariances, i.e., Cov/mu1/mu2:
nsig <- function(dist, f){
	N <- length(dist)-1
	
	#Some expressions occur repeatedly:
	mu <- mu.f(1, dist, f)
	mu10 <- mu.f(c(1,0), dist, f)
	nmu01mu11 <- N*mu.f(c(0,1), dist, f) - mu.f(c(1,1), dist, f)

	sig11 <- mu.f(2, dist, f)/mu^2 - 1
	sig12 <- mu.f(c(2,0), dist, f)/mu/mu10 - 1
	sig22 <- mu.f(c(2,0,0), dist, f)/mu10^2 - 1
	
	#The following differ from the unbounded case:
	sig13 <- (N*mu.f(c(1,1), dist, f) - mu.f(c(2,1), dist, f)) / mu / nmu01mu11 - 1
	sig23 <- (N*mu.f(c(1,0,1), dist, f) - mu.f(c(2,0,1), dist, f)) / mu10 / nmu01mu11 - 1
	sig33 <- (N^2*mu.f(c(0,1,1), dist, f) - 2*N*mu.f(c(1,1,1), dist, f) + mu.f(c(2,1,1), dist, f)) / nmu01mu11^2 - 1
	
	c(sig11, sig12, sig13, sig22, sig23, sig33)
}

#Asymptotic mean and SD:
TBSmeansd <- function(n, dist, f){
	N <- length(dist)-1
	mu <- mu.f(1, dist, f)
	
	Tf <- TBS(dist, f)
	nsigs <- nsig(dist, f)
	
	muTf <- Tf + Tf/n * sum(nsigs * c(N/(N-mu), -N/(N-mu), N/(N-mu), 0,-1,1))
	sigTf <- sqrt(Tf^2/n * sum(nsigs * c(N^2/(N-mu)^2, -2*N/(N-mu), 2*N/(N-mu), 1,-2,1)))
	
	c(muTf, sigTf)
}
TBSmeansd(156, dbinom(0:15, 15,0.15), function(x) exp(-x)) #1.0107784 0.1041094
TBSmeansd(156, dzib(0:15, 15,0.15, 0.4), function(x) exp(-x)) #0.43243590 0.05012159





#Testing H0: Bin with mu
#against H1: another dist

#Two-sided critical values:
crit.bin <- function(N,mu, n, f, level=0.05){
	z2 <- qnorm(1-level/2)
	pmf <- dbinom(0:N, N,mu/N)
	
	musig <- TBSmeansd(n, pmf, f)
	musig[1] + c(-1,1) * z2 * musig[2]
}
crit.bin(15,0.15, 156, function(x) exp(-x)) #0.8649775 1.1431412

#Exceedance probability for given interval "crit" under alternative:
apower <- function(n, dist, f, crit){
	N <- length(dist)-1
	
	musig <- TBSmeansd(n, dist, f)

	pnorm(crit[1], musig[1], musig[2]) + 1-pnorm(crit[2], musig[1], musig[2])
}
apower(156, dzib(0:15, 15,0.15, 0.02), function(x) exp(-x), crit.bin(15,0.15, 156, function(x) exp(-x))) #0.2361593





#############
#Power curves
#############




#Power with f(x)=u^x against ZIB alternative
#note that mu=(1-om)*N*p and Ibin=1+(N-1)*om*p/(1-(1-om)*p),
#note that mu=(1-om)*N*p and Ibin=1+(N-1)*om/(1-om)*mu/(N-mu).
#So N-Ibin = (N-1)*(1-p)/(1-mu/N)

N <- 15
mu <- 2
tabI <- c(1.1, 1.25, 1.5, 2)
tab.p <- 1-(N-tabI)/(N-1)*(1-mu/N)
tab.om <- 1-mu/N/tab.p
no.I <- length(tabI)
n <- 100

tab.u <- seq(0.01,0.99,0.01)
no.u <- length(tab.u)

tab.power <- array(NA, c(no.u, no.I))
for(k in 1:no.u){
	u <- tab.u[k]
	f <- function(x) u^x
	crit <- crit.bin(N,mu, n, f)
	
for(l in 1:no.I){
	p <- tab.p[l]
	om <- tab.om[l]
	pmf <- dzib(0:N, N,p, om)
	tab.power[k,l] <- apower(n, pmf, f, crit)
}}

cols <- grey(c(0,0.3,0.6))

matplot(tab.u, tab.power, type="l", lwd=2, ylim=c(0,1), col=cols, lty=1, xlab="u", ylab="Power", main=paste("H0: Bin; H1: ZIB; N=", N, ", mu=", mu, ", n=", n, sep=""))
legend("left",legend=tabI, title=expression(I[Bin]*"="), lwd=2, col=cols, lty=1, cex=0.8, bg=grey(1))
abline(h=0.05, lty=2)


#Find optimal u numerically:
for(l in 1:no.I){
	p <- tab.p[l]
	om <- tab.om[l]
	pmf <- dzib(0:N, N,p, om)

	f.optim <- function(u){
		f <- function(x) u^x
		-apower(n, pmf, f, crit.bin(N,mu, n, f))
	}
	
	u.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(0.9999), control=list(ndeps=c(1e-4)))
	print(c(n,mu,tabI[l], u.opt$par, -u.opt$value))
} #for I
#100   2.0000000   1.1000000   0.5880058   0.1559462
#100   2.0000000   1.2500000   0.4747804   0.6044684
#100   2.0000000   1.5000000   0.3048769   0.9946274






############
#Simulations
############

set.seed(123)
reps <- 1e5

N <- 10
mu <- 2
tabI <- c(1.1, 1.25, 1.5, 2)
no.I <- length(tabI)
tab.p <- 1-(N-tabI)/(N-1)*(1-mu/N)
tab.om <- 1-mu/N/tab.p
tabn <- c(50,100,250)
no.n <- length(tabn)



#Finite-sample performance of sizes:
reject <- array(0, c(no.I, no.n, 2))
opt.des <- array(0, c(no.I, no.n, 2))
for(l in 1:no.I){
	p <- tab.p[l]
	om <- tab.om[l]
	pmf <- dzib(0:N, N,p, om)

for(k in 1:no.n){
	n <- tabn[k]
	f.optim <- function(u){
		f <- function(x) u^x
		-apower(n, pmf, f, crit.bin(N,mu, n, f))
	}
	
	u.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(0.9999), control=list(ndeps=c(1e-4)))
	apow <- -u.opt$value
	u.opt <- u.opt$par
	opt.des[l,k,1] <- u.opt
	opt.des[l,k,2] <- apow

	f.opt <- function(x) u.opt^x
	crit.opt <- crit.bin(N,mu, n, f.opt)
	
	for(r in 1:reps){
		data0 <- rbinom(n, N,mu/N)
		data1 <- rzib(n, N,p, om)
		
		T0 <- sTBS(data0, N, f.opt)
		T1 <- sTBS(data1, N, f.opt)
		
		if(T0<crit.opt[1] || T0>crit.opt[2]) reject[l,k,1] <- reject[l,k,1]+1
		if(T1<crit.opt[1] || T1>crit.opt[2]) reject[l,k,2] <- reject[l,k,2]+1
	} #for reps
}} #for n,I

results <- array(NA, c(no.I*no.n, 8))
i <- 0
for(l in 1:no.I){
for(k in 1:no.n){
	i <- i+1
	results[i,] <- c(N,mu, tabI[l], tabn[k], opt.des[l,k,], reject[l,k,])
}} #for n,I
results

 # [1,]   10    2 1.10   50 0.7286171 0.1054077 4787  10641
 # [2,]   10    2 1.10  100 0.5910067 0.1565641 4914  15828
 # [3,]   10    2 1.10  250 0.4579141 0.3259475 5112  33216
 # [4,]   10    2 1.25   50 0.6108741 0.3428379 4908  34948
 # [5,]   10    2 1.25  100 0.4730473 0.6072824 5098  61693
 # [6,]   10    2 1.25  250 0.3043701 0.9628909 5173  95384
 # [7,]   10    2 1.50   50 0.4460435 0.8414264 5348  83273
 # [8,]   10    2 1.50  100 0.2971049 0.9953918 5555  98633
 # [9,]   10    2 1.50  250 0.3404277 1.0000000 5160 100000
# [10,]   10    2 2.00   50 0.2818484 0.9999985 6004  99924
# [11,]   10    2 2.00  100 0.6141239 1.0000000 4814 100000
# [12,]   10    2 2.00  250 0.9999000 1.0000000 4855 100000

 # [1,]   15    2 1.10   50 0.7205216 0.1050449 4853  10505
 # [2,]   15    2 1.10  100 0.5880058 0.1559462 4992  15808
 # [3,]   15    2 1.10  250 0.4595729 0.3240714 5100  33060
 # [4,]   15    2 1.25   50 0.6063961 0.3415751 4908  34771
 # [5,]   15    2 1.25  100 0.4747804 0.6044684 5058  61373
 # [6,]   15    2 1.25  250 0.3127851 0.9606893 5173  95240
 # [7,]   15    2 1.50   50 0.4501995 0.8370772 5339  83003
 # [8,]   15    2 1.50  100 0.3048769 0.9946274 5513  98571
 # [9,]   15    2 1.50  250 0.3082829 1.0000000 5238 100000
# [10,]   15    2 2.00   50 0.2883982 0.9999966 5923  99907
# [11,]   15    2 2.00  100 0.5928845 1.0000000 4868 100000
# [12,]   15    2 2.00  250 0.9999000 1.0000000 4917 100000



