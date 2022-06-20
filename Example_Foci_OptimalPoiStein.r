#Analyzing unbounded gamma-H2AX foci counts of Lopez et al. (2022)
#as discussed in Section 6.1.

#Absolute frequencies from zero to maximal count:

#Table 1. Frequency of foci observed for dose and post irradiation time (h).
#Dose=0, time in {0.5, 1, 2}
Time05 <- c(170,147,89,50,30,5,9)
Time1 <- c(286,81,67,52,12,2,0)
Time2 <- c(267,115,68,50,0,0,0)
data <- list(Time05, Time1, Time2)
no.data <- length(data)


#To be able to use existing costs, generate raw data first:
for(i in 1:no.data){
	x <- c()
	xmax <- length(data[[i]])-1
	for(k in 0:xmax) x <- c(x, rep(k, data[[i]][k+1]))
	data[[i]] <- x
	
	#The remaining code provides some descriptive statistics:
	n <- length(x)
	mu <- mean(x)	
	p0 <- mean(x==0)
	ID <- var(x)/mu
	
	U <- sqrt( (n-1)/2/(1-1/sum(x)) ) * (ID-1)

	print(c( n, mu, p0, ID, U ))
}
#500   1.348000   0.340000   1.447179   7.068696
#500   0.858000   0.572000   1.613759   9.705991
#500   0.802000   0.534000   1.287858   4.552553






##########################
#Stein tests for foci data
##########################

#Sample Stein-Chen statistic (3), data x, function f:
sTSC <- function(x, f){
	mean(x*f(x)) / mean(x) / mean(f(x+1))
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

#Population Stein-Chen statistic (4),
#for pmf vector "dist" and function f:
TSC <- function(dist, f){
	mu.f(c(1,0), dist, f) / mu.f(1, dist, f) / mu.f(c(0,1), dist, f)
}

#Asymptotics of test statistic according to Theorem 1:

#"Normalized" covariances, i.e., Cov/mu1/mu2,
#for pmf vector "dist" and function f:
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

#Asymptotic mean and SD with sample size n,
#for pmf vector "dist" and function f:
TSCmeansd <- function(n, dist, f){
	Tf <- TSC(dist, f)
	nsigs <- nsig(dist, f)
	
	muTf <- Tf + Tf/n * sum(nsigs * c(1,-1,1,0,-1,1))
	sigTf <- sqrt(Tf^2/n * sum(nsigs * c(1,-2,2,1,-2,1)))
	
	c(muTf, sigTf)
}


#Testing H0: Poi with mean mu
#against H1: another dist

#Two-sided critical values for given level,
#with mean mu, sample size n, function f, and upper truncation:
crit.poi <- function(mu, n, f, upper=50, level=0.05){
	z2 <- qnorm(1-level/2)
	pmf <- dpois(0:upper, mu)
	
	musig <- TSCmeansd(n, pmf, f)
	musig[1] + c(-1,1) * z2 * musig[2]
}

#Exceedance probability for given interval "crit" under alternative (specified by pmf vector "dist"),
#with sample size n and function f,
#where crit might be computed using function "crit.poi":
apower <- function(n, dist, f, crit){
	musig <- TSCmeansd(n, dist, f)

	pnorm(crit[1], musig[1], musig[2]) + 1-pnorm(crit[2], musig[1], musig[2])
}

#Asymptotic p value of statistic TSC,
#with mean mu, sample size n, function f, and upper truncation:
apvalue <- function(TSC, mu, n, f, upper=50){
	pmf <- dpois(0:upper, mu)
	
	musig <- TSCmeansd(n, pmf, f)

	2*(1-pnorm(abs(TSC-musig[1])/musig[2]))
}
#apvalue(TSC, mu, n, f)



#####################################
#data tested against NB alternative
#####################################


#We use individual means of data for test design.
#The following code generates Table 8 (except competitor tests):

upper <- 50

tabI <- c(1.1, 1.25, 1.5)
no.I <- length(tabI)

#Find optimal b numerically:
for(i in 1:no.data){
	x <- data[[i]]
	mu <- mean(x)
	n <- length(x)

	tabnu <- mu/(tabI-1)

for(l in 1:no.I){
	nu <- tabnu[l]
	pmf <- dnbinom(0:upper, nu, nu/(nu+mu))

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
} #for data
#500   1.3480000   1.1000000   1.0579878   0.3693143   1.3369225   0.9033330   1.0934839   0.0000000
#500   1.3480000   1.2500000   0.8731352   0.9288508   1.3128871   0.9129810   1.0844950   0.0000000
#500   1.3480000   1.5000000   0.6955811   0.9997889   1.2882158   0.9203027   1.0777408   0.0000000
#500   0.8580000   1.1000000   1.0945090   0.3721618   1.6909806   0.8466755   1.1480820   0.0000000
#500   0.8580000   1.2500000   0.8481146   0.9205671   1.7438653   0.8619597   1.1342340   0.0000000
#500   0.8580000   1.5000000   0.6280213   0.9995964   1.7884310   0.8707548   1.1265484   0.0000000
#500   0.8020000   1.1000000   1.1014859   0.3727153   1.3320052   0.8357312   1.1586095   0.0000480
#500   0.8020000   1.2500000   0.8431795   0.9190112   1.3905248   0.8521198   1.1438405   0.0000001
#500   0.8020000   1.5000000   0.6162887   0.9995493   1.4364632   0.8610518   1.1361313   0.0000000







####################
#Further Simulations
#for Table 8
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

tabmu <- tabn <- rep(NA, no.data)
for(i in 1:no.data){
	x <- data[[i]]
	tabmu[i] <- mean(x)
	tabn[i] <- length(x)
}
tabI <- c(1.1, 1.25, 1.5)
no.I <- length(tabI)



#Finite-sample performance of sizes:
reject <- array(0, c(no.I, no.data, 2))
for(k in 1:no.data){
	n <- tabn[k]
	mu <- tabmu[k]
	
for(l in 1:no.I){
	Id <- tabI[l]
	nu <- mu/(Id-1)


	for(r in 1:reps){
		data0 <- rpois(n, mu)
		data1 <- rnbinom(n, nu, nu/(nu+mu))
		
		T0 <- Upoi(data0)
		T1 <- Upoi(data1)
		
		if(T0<crit[1] || T0>crit[2]) reject[l,k,1] <- reject[l,k,1]+1
		if(T1<crit[1] || T1>crit[2]) reject[l,k,2] <- reject[l,k,2]+1
	} #for reps

	#print(c(mu,Id,n, opt.des[l,k,], reject[l,k,]))
}} #for no.data,I

results <- array(NA, c(no.I*no.data, 5))
i <- 0
for(k in 1:no.data){
for(l in 1:no.I){
	i <- i+1
	results[i,] <- c(tabmu[k], tabI[l], tabn[k], reject[l,k,])
}} #for no.data,I
results
 # [1,] 1.348 1.10  500 4942 35785
 # [2,] 1.348 1.25  500 4915 93348
 # [3,] 1.348 1.50  500 4966 99996
 # [4,] 0.858 1.10  500 4970 36147
 # [5,] 0.858 1.25  500 4999 92554
 # [6,] 0.858 1.50  500 5013 99989
 # [7,] 0.802 1.10  500 4961 35539
 # [8,] 0.802 1.25  500 5020 92263
 # [9,] 0.802 1.50  500 4883 99987

#Power values of last column used in Table 8.





#For comparison, try out discrete KS-test, 
#as implemented in dgof package, see Arnold and Emerson (2011).
#Null hypothesis needs to be specified as stepfun.
#Exact P-values not offered for sample size n>30.

#install.packages("dgof")
library("dgof")

level <- 0.05 #actual target level


set.seed(123)
reps <- 1e5

upper <- 50

tabmu <- tabn <- rep(NA, no.data)
for(i in 1:no.data){
	x <- data[[i]]
	tabmu[i] <- mean(x)
	tabn[i] <- length(x)
}
tabI <- c(1.1, 1.25, 1.5)
no.I <- length(tabI)



#Finite-sample performance of sizes:
reject <- array(0, c(no.I, no.data, 2))
for(k in 1:no.data){
	n <- tabn[k]
	mu <- tabmu[k]

	#CDF under null:
	f0 <- stepfun(0:upper, ppois((-1):upper, mu))
	
for(l in 1:no.I){
	Id <- tabI[l]
	nu <- mu/(Id-1)


	for(r in 1:reps){
		data0 <- rpois(n, mu)
		data1 <- rnbinom(n, nu, nu/(nu+mu))
		
		P0 <- ks.test(data0, f0)$p.value
		P1 <- ks.test(data1, f0)$p.value
		
		if(P0<level) reject[l,k,1] <- reject[l,k,1]+1
		if(P1<level) reject[l,k,2] <- reject[l,k,2]+1
	} #for reps

	#print(c(mu,Id,n, opt.des[l,k,], reject[l,k,]))
}} #for no.data,I

results <- array(NA, c(no.I*no.data, 5))
i <- 0
for(k in 1:no.data){
for(l in 1:no.I){
	i <- i+1
	results[i,] <- c(tabmu[k], tabI[l], tabn[k], reject[l,k,])
}} #for no.data,I
results
 # [1,] 1.348 1.10  500  712  2065
 # [2,] 1.348 1.25  500  738 16688
 # [3,] 1.348 1.50  500  729 76350
 # [4,] 0.858 1.10  500  635  2620
 # [5,] 0.858 1.25  500  632 18546
 # [6,] 0.858 1.50  500  617 73776
 # [7,] 0.802 1.10  500  627  2710
 # [8,] 0.802 1.25  500  647 18420
 # [9,] 0.802 1.50  500  669 72155

#These results are summarized in Table S9 of the Supplement.




#For comparison, try out discrete KS-test, 
#as implemented in dgof package, see Arnold and Emerson (2011),
#but where critical value determined by simulation.

#install.packages("dgof")
library("dgof")

level <- 0.05 #actual target level


set.seed(123)
reps <- 1e5

upper <- 50

tabmu <- tabn <- rep(NA, no.data)
for(i in 1:no.data){
	x <- data[[i]]
	tabmu[i] <- mean(x)
	tabn[i] <- length(x)
}
tabI <- c(1.1, 1.25, 1.5)
no.I <- length(tabI)



#Finite-sample performance of sizes:
reject <- array(0, c(no.I, no.data, 2))
tabcrit <- rep(NA, no.data)
for(k in 1:no.data){
	n <- tabn[k]
	mu <- tabmu[k]

	#CDF under null:
	f0 <- stepfun(0:upper, ppois((-1):upper, mu))

	#Determine critical value by simulation:
	D <- rep(NA, reps)
	for(r in 1:reps){
		data0 <- rpois(n, mu)
		D[r] <- ks.test(data0, f0)$statistic
	} #for reps
	tabcrit[k] <- quantile(D, 1-level)[[1]]
	
for(l in 1:no.I){
	Id <- tabI[l]
	nu <- mu/(Id-1)


	for(r in 1:reps){
		data0 <- rpois(n, mu)
		data1 <- rnbinom(n, nu, nu/(nu+mu))
		
		T0 <- ks.test(data0, f0)$statistic
		T1 <- ks.test(data1, f0)$statistic
		
		if(T0>tabcrit[k]) reject[l,k,1] <- reject[l,k,1]+1
		if(T1>tabcrit[k]) reject[l,k,2] <- reject[l,k,2]+1
	} #for reps
}} #for no.data,I

results <- array(NA, c(no.I*no.data, 6))
i <- 0
for(k in 1:no.data){
for(l in 1:no.I){
	i <- i+1
	results[i,] <- c(tabmu[k], tabI[l], tabn[k], tabcrit[k], reject[l,k,])
}} #for no.data,I
results
 # [1,] 1.348 1.10  500 0.04608526 5015 11379
 # [2,] 1.348 1.25  500 0.04608526 4805 45003
 # [3,] 1.348 1.50  500 0.04608526 4972 95091
 # [4,] 0.858 1.10  500 0.04587937 5078 12158
 # [5,] 0.858 1.25  500 0.04587937 5169 44588
 # [6,] 0.858 1.50  500 0.04587937 5181 93168
 # [7,] 0.802 1.10  500 0.04556880 4716 10855
 # [8,] 0.802 1.25  500 0.04556880 4711 41360
 # [9,] 0.802 1.50  500 0.04556880 4764 91475

#These results are summarized in Table S9 of the Supplement.

#Power values of last column used in Table 8.





#For comparison, try out LR-test, 
#with asymptotic design.

#Negative log-likelihood function of NB distribution for data x,
#if using parametrization by mean and dispersion index:
llNB <- function(par, x){
	#par = c(mu,Id)
	mu <- par[1]
	Id <- par[2]
	-sum(log(dnbinom(x, mu/(Id-1), 1/Id)))
}

#Function to execute the LR-test for Poi-H0 against NB-H1:
LRpoinb <- function(x){
	#MLE of Poisson mean:
	mu.poi <- mean(x)
	lmax.poi <- sum(log(dpois(x, mu.poi)))
	
	#MLE of NB distribution by numerical optimization:
	estml <- suppressWarnings(optim(c(mu.poi,max(1.1,var(x)/mu.poi)), llNB, method="L-BFGS-B", lower=c(0.001,1.001), upper=c(9.999,9.999), control=list(ndeps=c(1e-4,1e-4)), x=x, hessian=FALSE))
	
	mu.nb <- estml$par[1]
	Id.nb <- estml$par[2]
	lmax.nb <- (-estml$value)

	c(mu.poi, lmax.poi, mu.nb, Id.nb, lmax.nb, 2*(lmax.nb-lmax.poi))
	#The last value in this output is the actual test statistic.
}


level <- 0.05 #actual target level
crit <- qchisq(1-2*level, 1) #(Self & Liang, 1987, Case 5)

set.seed(123)
reps <- 1e5

tabmu <- tabn <- rep(NA, no.data)
for(i in 1:no.data){
	x <- data[[i]]
	tabmu[i] <- mean(x)
	tabn[i] <- length(x)
}
tabI <- c(1.1, 1.25, 1.5)
no.I <- length(tabI)



#Finite-sample performance of sizes:
reject <- array(0, c(no.I, no.data, 2))
for(k in 1:no.data){
	n <- tabn[k]
	mu <- tabmu[k]
	
for(l in 1:no.I){
	Id <- tabI[l]
	nu <- mu/(Id-1)


	for(r in 1:reps){
		data0 <- rpois(n, mu)
		data1 <- rnbinom(n, nu, nu/(nu+mu))
		
		T0 <- LRpoinb(data0)[6]
		T1 <- LRpoinb(data1)[6]
		
		if(T0>crit) reject[l,k,1] <- reject[l,k,1]+1
		if(T1>crit) reject[l,k,2] <- reject[l,k,2]+1
	} #for reps

	#print(c(mu,Id,n, opt.des[l,k,], reject[l,k,]))
}} #for no.data,I

results <- array(NA, c(no.I*no.data, 5))
i <- 0
for(k in 1:no.data){
for(l in 1:no.I){
	i <- i+1
	results[i,] <- c(tabmu[k], tabI[l], tabn[k], reject[l,k,])
}} #for no.data,I
results
# 1.348 1.1 500 4450 42643
# 1.348 1.25 500 4390 95355
# 1.348 1.5 500 4558 99999
# 0.858 1.1 500 4342 42273
# 0.858 1.25 500 4512 94581
# 0.858 1.5 500 4326 99994
# 0.802 1.1 500 4436 41803
# 0.802 1.25 500 4511 94401
# 0.802 1.5 500 4269 99995

#Power values of last column used in Table 8.


