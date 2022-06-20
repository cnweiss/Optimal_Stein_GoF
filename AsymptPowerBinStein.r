#Stein statistic for Bin-null.

#Simulation results of Section 5 and Supplement S1,
#where asymptotic power computed by using Theorem 2 in Section 4.


#Functions for ZIB distribution with
#default parametrization (p,omega) and upper bound N:

#Probability mass function in x:
dzib <- function(x, N,p, omega=0) omega*(x==0) + (1-omega)*dbinom(x, N,p)

#Cumulative distribution function in x:
pzib <- function(x, N,p, omega=0) omega*(x>=0) + (1-omega)*pbinom(x, N,p)

#n i.i.d. random numbers:
rzib <- function(n, N,p, omega=0) rbinom(n, 1, 1-omega)*rbinom(n, N,p)

#Quantile function in quantile level ql:
qzib <- function(ql, N,p, omega=0){
	up <- qbinom(ql, N,p) #ZIB quantiles not later than Bin quantiles
	cdf <- pzib(0:up, N,p, omega)
	min((0:up)[cdf>=ql])
}
qzib <- Vectorize(qzib)
#These functions are used for the alternative hypothesis.



#Now, Stein statistics and asymptotics from manuscript.

#Sample Bin-Stein statistic (14), data x, upper bound N, function f:
sTBS <- function(x, N, f){
	(N/mean(x)-1) * mean(x*f(x)) / mean((N-x)*f(x+1))
}


#Exact moment computation (16) of order "klvec" and function f, 
#where upper bound N implicitly provided by length of pmf vector "dist":
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
#mu.f(c(0,1), dbinom(0:15, 15,0.15), function(x) exp(-x)) #0.08255503
#mu.f(c(2), dbinom(0:15, 15,0.15), function(x) exp(-x)) #6.975

#Population Bin-Stein statistic (15),
#for pmf vector "dist" and function f:
TBS <- function(dist, f){
	N <- length(dist)-1

	(N/mu.f(1, dist, f)-1) * mu.f(c(1,0), dist, f) / (N*mu.f(c(0,1), dist, f) - mu.f(c(1,1), dist, f))
}
#TBS(dbinom(0:15, 15,0.15), function(x) exp(-x)) #1
#TBS(dzib(0:15, 15,0.15, 0.4), function(x) exp(-x)) #0.4285461



#Asymptotics of test statistics according to Theorem 2:

#"Normalized" covariances, i.e., Cov/mu1/mu2,
#for pmf vector "dist" and function f:
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

#Asymptotic mean and SD with sample size n,
#for pmf vector "dist" and function f:
TBSmeansd <- function(n, dist, f){
	N <- length(dist)-1
	mu <- mu.f(1, dist, f)
	
	Tf <- TBS(dist, f)
	nsigs <- nsig(dist, f)
	
	muTf <- Tf + Tf/n * sum(nsigs * c(N/(N-mu), -N/(N-mu), N/(N-mu), 0,-1,1))
	sigTf <- sqrt(Tf^2/n * sum(nsigs * c(N^2/(N-mu)^2, -2*N/(N-mu), 2*N/(N-mu), 1,-2,1)))
	
	c(muTf, sigTf)
}
#TBSmeansd(156, dbinom(0:15, 15,0.15), function(x) exp(-x)) #1.0107784 0.1041094
#TBSmeansd(156, dzib(0:15, 15,0.15, 0.4), function(x) exp(-x)) #0.43243590 0.05012159





#Testing H0: Bin with mean mu
#against H1: another dist,
#where upper bound N


#Two-sided critical values for given level,
#with mean mu, sample size n, function f, and upper bound N:
crit.bin <- function(N,mu, n, f, level=0.05){
	z2 <- qnorm(1-level/2)
	pmf <- dbinom(0:N, N,mu/N)
	
	musig <- TBSmeansd(n, pmf, f)
	musig[1] + c(-1,1) * z2 * musig[2]
}
#crit.bin(15,0.15, 156, function(x) exp(-x)) #0.8649775 1.1431412

#Exceedance probability for given interval "crit" under alternative (specified by pmf vector "dist"),
#with sample size n and function f,
#where crit might be computed using function "crit.bin":
apower <- function(n, dist, f, crit){
	N <- length(dist)-1
	
	musig <- TBSmeansd(n, dist, f)

	pnorm(crit[1], musig[1], musig[2]) + 1-pnorm(crit[2], musig[1], musig[2])
}
#apower(156, dzib(0:15, 15,0.15, 0.02), function(x) exp(-x), crit.bin(15,0.15, 156, function(x) exp(-x))) #0.2361593





######################
#Example code for 
#creating power curves
#(not in manuscript)
######################


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






###########################
#Simulations for Section 5:
###########################

#ZIB alternative with f(x)=|x-1|^a, see Table 5:

set.seed(123)
reps <- 1e5

N <- 10 #15, switch between both choices of N
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
	f.optim <- function(a){
		f <- function(x) abs(x-1)^a
		-apower(n, pmf, f, crit.bin(N,mu, n, f))
	}
	
	a.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(4.9999), control=list(ndeps=c(1e-4)))
	apow <- -a.opt$value
	a.opt <- a.opt$par
	opt.des[l,k,1] <- a.opt
	opt.des[l,k,2] <- apow

	f.opt <- function(x) abs(x-1)^a.opt
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
 # [1,]   10    2 1.10   50 0.2361522 0.1593182 5493  15433
 # [2,]   10    2 1.10  100 0.2328368 0.2275723 5195  21759
 # [3,]   10    2 1.10  250 0.2296558 0.4157375 5127  40052
 # [4,]   10    2 1.25   50 0.2252030 0.5084546 5465  46575
 # [5,]   10    2 1.25  100 0.2228796 0.7233589 5225  70522
 # [6,]   10    2 1.25  250 0.2198441 0.9594404 5076  96725
 # [7,]   10    2 1.50   50 0.2103355 0.8814632 5592  88867
 # [8,]   10    2 1.50  100 0.2078270 0.9814048 5387  99140
 # [9,]   10    2 1.50  250 0.2041088 0.9999502 5106 100000
# [10,]   10    2 2.00   50 0.1827318 0.9917711 5393  99923
# [11,]   10    2 2.00  100 0.1794547 0.9998505 5269 100000
# [12,]   10    2 2.00  250 0.9999971 0.9999997 4890 100000

 # [1,]   15    2 1.10   50 0.2393027 0.1547187 5565  15155
 # [2,]   15    2 1.10  100 0.2358677 0.2217736 5224  21346
 # [3,]   15    2 1.10  250 0.2326149 0.4078278 5144  39418
 # [4,]   15    2 1.25   50 0.2283645 0.4976242 5496  45533
 # [5,]   15    2 1.25  100 0.2259819 0.7149003 5286  69654
 # [6,]   15    2 1.25  250 0.2228990 0.9575011 5059  96557
 # [7,]   15    2 1.50   50 0.2138365 0.8763131 5608  88169
 # [8,]   15    2 1.50  100 0.2112650 0.9805091 5288  99064
 # [9,]   15    2 1.50  250 0.2075011 0.9999487 5150 100000
# [10,]   15    2 2.00   50 0.1875531 0.9913005 5558  99909
# [11,]   15    2 2.00  100 0.1842129 0.9998447 5200 100000
# [12,]   15    2 2.00  250 0.9999957 0.9999996 4966 100000





#ZIB alternative with f(x)=u^x, see Table 6:

set.seed(123)
reps <- 1e5

N <- 10 #15, switch between both choices of N
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









###############################
#Simulations for Supplement S1:
###############################

#ZIB alternative with f(x)=(x^b-1)/b, see Table S6:

set.seed(123)
reps <- 1e5

N <- 10 #15, switch between both choices of N
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
	f.optim <- function(b){
		f <- function(x) (x^b-1)/b
		-apower(n, pmf, f, crit.bin(N,mu, n, f))
	}
	
	b.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(4.9999), control=list(ndeps=c(1e-4)))
	apow <- -b.opt$value
	b.opt <- b.opt$par
	opt.des[l,k,1] <- b.opt
	opt.des[l,k,2] <- apow

	f.opt <- function(x) (x^b.opt-1)/b.opt
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
 # [1,]   10    2 1.10   50 0.0001000 0.1399452 5330  13679
 # [2,]   10    2 1.10  100 0.0001000 0.1950810 5129  18693
 # [3,]   10    2 1.10  250 0.0001000 0.3519384 5077  33917
 # [4,]   10    2 1.25   50 0.0001000 0.4441222 5307  40329
 # [5,]   10    2 1.25  100 0.0001000 0.6498184 5123  62505
 # [6,]   10    2 1.25  250 0.0001000 0.9256932 5053  93235
 # [7,]   10    2 1.50   50 0.0001000 0.8366688 5413  83266
 # [8,]   10    2 1.50  100 0.0001000 0.9661522 5314  97837
 # [9,]   10    2 1.50  250 0.0001000 0.9997869 5091 100000
# [10,]   10    2 2.00   50 0.0001000 0.9864950 5181  99786
# [11,]   10    2 2.00  100 0.0001000 0.9996376 4988 100000
# [12,]   10    2 2.00  250 0.9999985 0.9999997 4890 100000

 # [1,]   15    2 1.10   50 0.0001000 0.1359991 5301  13358
 # [2,]   15    2 1.10  100 0.0001000 0.1900574 5190  18420
 # [3,]   15    2 1.10  250 0.0001000 0.3446515 5094  33350
 # [4,]   15    2 1.25   50 0.0001000 0.4330675 5273  39499
 # [5,]   15    2 1.25  100 0.0001000 0.6397981 5162  61515
 # [6,]   15    2 1.25  250 0.0001000 0.9220410 4994  92906
 # [7,]   15    2 1.50   50 0.0001000 0.8291970 5301  82278
 # [8,]   15    2 1.50  100 0.0001000 0.9642118 5285  97640
 # [9,]   15    2 1.50  250 0.0001000 0.9997730 5129  99999
# [10,]   15    2 2.00   50 0.0001000 0.9855568 5311  99750
# [11,]   15    2 2.00  100 0.0001000 0.9996141 4956 100000
# [12,]   15    2 2.00  250 0.9999978 0.9999996 4966 100000





#ZIB alternative with f(x)=1/(x+1)^d, see Table S7:

set.seed(123)
reps <- 1e5

N <- 10 #15, switch between both choices of N
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
	f.optim <- function(d){
		f <- function(x) 1/(x+1)^d
		-apower(n, pmf, f, crit.bin(N,mu, n, f))
	}
	
	d.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(3), control=list(ndeps=c(1e-4)))
	apow <- -d.opt$value
	d.opt <- d.opt$par
	opt.des[l,k,1] <- d.opt
	opt.des[l,k,2] <- apow

	f.opt <- function(x) 1/(x+1)^d.opt
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
 # [1,]   10    2 1.10   50 0.3853052 0.1153097 4761  11477
 # [2,]   10    2 1.10  100 0.8878371 0.1700602 4842  17007
 # [3,]   10    2 1.10  250 1.5258834 0.3437381 5066  34604
 # [4,]   10    2 1.25   50 0.7938694 0.3764186 4811  37340
 # [5,]   10    2 1.25  100 1.4749850 0.6338039 4988  63796
 # [6,]   10    2 1.25  250 2.7242265 0.9642995 5086  95695
 # [7,]   10    2 1.50   50 1.7590406 0.8556986 5116  84844
 # [8,]   10    2 1.50  100 2.9658139 0.9954459 5372  98727
 # [9,]   10    2 1.50  250 1.0000087 0.9999979 5014 100000
# [10,]   10    2 2.00   50 3.0000000 0.9999982 5633  99935
# [11,]   10    2 2.00  100 1.0000001 1.0000000 4793 100000
# [12,]   10    2 2.00  250 1.0000000 1.0000000 4891 100000

 # [1,]   15    2 1.10   50 0.4377007 0.1142201 4819  11379
 # [2,]   15    2 1.10  100 0.9233177 0.1686206 4910  16930
 # [3,]   15    2 1.10  250 1.5366033 0.3409782 5027  34370
 # [4,]   15    2 1.25   50 0.8375353 0.3731620 4788  37112
 # [5,]   15    2 1.25  100 1.4863149 0.6298172 4943  63358
 # [6,]   15    2 1.25  250 2.6610619 0.9621868 5070  95577
 # [7,]   15    2 1.50   50 1.7484483 0.8509367 5048  84396
 # [8,]   15    2 1.50  100 2.9041949 0.9946901 5355  98646
 # [9,]   15    2 1.50  250 1.0000114 0.9999972 5043 100000
# [10,]   15    2 2.00   50 3.0000000 0.9999961 5586  99915
# [11,]   15    2 2.00  100 1.0000001 1.0000000 4766 100000
# [12,]   15    2 2.00  250 1.0000000 1.0000000 5002 100000



