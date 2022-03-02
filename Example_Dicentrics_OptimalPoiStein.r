#Analyzing unbounded dicentrics counts of Pujol et al. (2016)

#Absolute frequencies from zero to maximal count:

#Table 1B. Heterogeneous simulations (mixed-Poisson alternative)
M4F1 <- c(275, 66, 42, 15, 0, 1)
F4M1 <- c(208, 61, 34, 17, 5, 1)
M6F2 <- c(166, 76, 23, 13, 10)
F6M2 <- c(160, 55, 19, 17, 9, 4)
dataB <- list(M4F1, F4M1, M6F2, F6M2)
nB <- length(dataB)


#Functions for mixed-Poi distribution

#Default parametrization (lambda1, lambda2, omega):

dmpoi <- function(x, lambda1, lambda2, omega=0.5) omega*dpois(x, lambda1) + (1-omega)*dpois(x, lambda2)

pmpoi <- function(x, lambda1, lambda2, omega=0.5) omega*ppois(x, lambda1) + (1-omega)*ppois(x, lambda2)

rmpoi <- function(n, lambda1, lambda2, omega=0.5){
	X <- rep(NA, n)
	W <- rbinom(n, 1, 1-omega)
	n2 <- sum(W)
	X[W==1] <- rpois(n2, lambda2)
	X[W==0] <- rpois(n-n2, lambda1)
	X
}






for(i in 1:nB){
	#Generate raw data first:
	x <- c()
	xmax <- length(dataB[[i]])-1
	for(k in 0:xmax) x <- c(x, rep(k, dataB[[i]][k+1]))
	dataB[[i]] <- x
	
	n <- length(x)
	mu <- mean(x)	
	p0 <- mean(x==0)
	ID <- var(x)/mu
	
	U <- sqrt( (n-1)/2/(1-1/sum(x)) ) * (ID-1)

	print(c( n, mu, p0, ID, U ))
}
#399   0.5012531   0.6892231   1.4724372   6.6812708
#326   0.6288344   0.6380368   1.5955722   7.6106714
#288   0.6979167   0.5763889   1.5212960   6.2602740
#264   0.7575758   0.6060606   1.8895817  10.2267397






################################
#Stein tests for dicentrics data
################################



#Computing asymptotic power curves using Theorem 4.1.



#Sample Stein-Chen statistic:
sTSC <- function(x, f){
	mean(x*f(x)) / mean(x) / mean(f(x+1))
}


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

#True Stein-Chen statistic:
TSC <- function(dist, f){
	mu.f(c(1,0), dist, f) / mu.f(1, dist, f) / mu.f(c(0,1), dist, f)
}

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


#Testing H0: Poi with mu
#against H1: another dist

#Two-sided critical values:
crit.poi <- function(mu, n, f, upper=50, level=0.05){
	z2 <- qnorm(1-level/2)
	pmf <- dpois(0:upper, mu)
	
	musig <- TSCmeansd(n, pmf, f)
	musig[1] + c(-1,1) * z2 * musig[2]
}

#Exceedance probability for given interval "crit" under alternative:
apower <- function(n, dist, f, crit){
	musig <- TSCmeansd(n, dist, f)

	pnorm(crit[1], musig[1], musig[2]) + 1-pnorm(crit[2], musig[1], musig[2])
}

#Asymptotic p value:
apvalue <- function(TSC, mu, n, f, upper=50){
	pmf <- dpois(0:upper, mu)
	
	musig <- TSCmeansd(n, pmf, f)

	2*(1-pnorm(abs(TSC-musig[1])/musig[2]))
}
#apvalue(TSC, mu, n, f)






###########################################
#dataB tested against mixed-Poi alternative
###########################################


#Use individual means of data for test design:

upper <- 50

tabI <- c(1.1, 1.25, 1.5)
no.I <- length(tabI)

tab.a <- seq(0.1,3,0.1)
no.a <- length(tab.a)

#Find optimal b numerically:
for(i in 1:nB){
	x <- dataB[[i]]
	mu <- mean(x)
	n <- length(x)

	tablam1 <- mu - sqrt((tabI-1)*mu)
	tablam2 <- mu + sqrt((tabI-1)*mu)

for(l in 1:no.I){
	lam1 <- tablam1[l]
	lam2 <- tablam2[l]
	pmf <- dmpoi(0:upper, lam1, lam2)

	f.optim <- function(a){
		f <- function(x) abs(x-1)^a
		-apower(n, pmf, f, crit.poi(mu, n, f))
	}
	
	a.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(4.9999), control=list(ndeps=c(1e-4)))
	apow <- -a.opt$value
	a.opt <- a.opt$par
	f.opt <- function(x) abs(x-1)^a.opt
	print(round(c(n,mu,tabI[l], a.opt, apow, sTSC(x, f.opt), crit.poi(mu, n, f.opt), apvalue(sTSC(x, f.opt), mu, n, f.opt)), 7))
} #for I
} #for dataB
#399   0.5012531   1.10   0.9511542   0.3191954   1.9459646   0.7221545   1.2684758   0.0000000
#399   0.5012531   1.25   0.5706331   0.8999999   2.0314194   0.7415645   1.2533563   0.0000000
#399   0.5012531   1.50   0.3186740   0.9999390   2.0876031   0.7447543   1.2525623   0.0000000
#326   0.6288344   1.10   0.9877132   0.2769839   1.9398007   0.7520818   1.2383169   0.0000000
#326   0.6288344   1.25   0.6438949   0.8436881   1.9481590   0.7737111   1.2205552   0.0000000
#326   0.6288344   1.50   0.3872297   0.9995944   1.9484111   0.7796099   1.2171135   0.0000000
#288   0.6979167   1.10   1.0061196   0.2546535   1.7405744   0.7604297   1.2295414   0.0000000
#288   0.6979167   1.25   0.6763451   0.8034489   1.6646262   0.7834143   1.2103779   0.0000000
#288   0.6979167   1.50   0.4239796   0.9988790   1.5950949   0.7908245   1.2054656   0.0000000
#264   0.7575758   1.10   1.0176232   0.2402484   2.1683038   0.7681293   1.2216432   0.0000000
#264   0.7575758   1.25   0.6996156   0.7732972   2.0909682   0.7917661   1.2017316   0.0000000
#264   0.7575758   1.50   0.4539894   0.9978579   2.0129376   0.8002665   1.1956986   0.0000000





####################
#Further Simulations
####################

#For comparison, try out Poisson U-test,
#see Savage (1970, p. 188), R&C (1956, p. 266), and
#IAEA (2011, p. 48) [EPR-Biodosimetry 2011_web.pdf],
#because U-test also used in Pujol.

Upoi <- function(x){
	sqrt( (length(x)-1)/2/(1-1/sum(x)) ) * (var(x)/mean(x)-1)
}
level <- 0.05
crit <- c(-1,1) * qnorm(1-level/2)


set.seed(123)
reps <- 1e5

upper <- 50

tabmu <- tabn <- rep(NA, nB)
for(i in 1:nB){
	x <- dataB[[i]]
	tabmu[i] <- mean(x)
	tabn[i] <- length(x)
}
tabI <- c(1.1, 1.25, 1.5)
no.I <- length(tabI)



#Finite-sample performance of sizes:
reject <- array(0, c(no.I, nB, 2))
for(k in 1:nB){
	n <- tabn[k]
	mu <- tabmu[k]
	
for(l in 1:no.I){
	Id <- tabI[l]
	lam1 <- mu - sqrt((Id-1)*mu)
	lam2 <- mu + sqrt((Id-1)*mu)

	for(r in 1:reps){
		data0 <- rpois(n, mu)
		data1 <- rmpoi(n, lam1, lam2)
		
		T0 <- Upoi(data0)
		T1 <- Upoi(data1)
		
		if(T0<crit[1] || T0>crit[2]) reject[l,k,1] <- reject[l,k,1]+1
		if(T1<crit[1] || T1>crit[2]) reject[l,k,2] <- reject[l,k,2]+1
	} #for reps

	#print(c(mu,Id,n, opt.des[l,k,], reject[l,k,]))
}} #for n,I

results <- array(NA, c(no.I*nB, 5))
i <- 0
for(k in 1:nB){
for(l in 1:no.I){
	i <- i+1
	results[i,] <- c(tabmu[k], tabI[l], tabn[k], reject[l,k,])
}} #for n,I
results
 # [1,] 0.5012531 1.10  399 4875 30514
 # [2,] 0.5012531 1.25  399 5088 88119
 # [3,] 0.5012531 1.50  399 4960 99995
 # [4,] 0.6288344 1.10  326 4919 26668
 # [5,] 0.6288344 1.25  326 4957 82172
 # [6,] 0.6288344 1.50  326 4933 99956
 # [7,] 0.6979167 1.10  288 4946 24349
 # [8,] 0.6979167 1.25  288 4921 77871
 # [9,] 0.6979167 1.50  288 4902 99884
# [10,] 0.7575758 1.10  264 4944 22967
# [11,] 0.7575758 1.25  264 4939 74547
# [12,] 0.7575758 1.50  264 4855 99782
