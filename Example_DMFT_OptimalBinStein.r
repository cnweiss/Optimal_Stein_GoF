#Analyzing bounded DMFT counts of Böhning et al. (1999)
# Variables:
# dmftb - DMFT-Index at the beginning of the study
# dmfte - DMFT-Index at the end of the study
# gender - gender of the school child
# 0 - female
# 1 - male
# ethnic - ethnic group
# 1 - dark
# 2 - white
# 3 - black
# school - school (=kind of prevention measure)
# 1 - oral health education
# 2 - all four methods together
# 3 - control school (no prevention measure)
# 4 - enrichment of the school diet with ricebran
# 5 - mouthrinse with 0.2% NaF-solution
# 6 - oral hygiene


data <- read.table(file=file.choose(), col.names=c("dmftb", "dmfte", "gender", "ethnic", "school")) #dmft_bohning99.txt
n <- dim(data)[1]
n #797



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


#negative log-likelihood for ZIB:
llzib <- function(par, x, N){
	#par = c(p,om)
	p <- par[1]
	om <- par[2]
	
	-sum(log(dzib(x, N,p,om)))
}

#eight deciduous molars considered, hence
N <- 8






###########################
#Stein tests for dmfte data
###########################



#Computing asymptotic power curves using Theorem 4.1.



#Sample Bin-Stein statistic:
sTBS <- function(x, N, f){
	(N/mean(x)-1) * mean(x*f(x)) / mean((N-x)*f(x+1))
}

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

#True Bin-Stein statistic:
TBS <- function(dist, f){
	N <- length(dist)-1

	(N/mu.f(1, dist, f)-1) * mu.f(c(1,0), dist, f) / (N*mu.f(c(0,1), dist, f) - mu.f(c(1,1), dist, f))
}

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


#Testing H0: Bin with mu
#against H1: another dist

#Two-sided critical values:
crit.bin <- function(N,mu, n, f, level=0.05){
	z2 <- qnorm(1-level/2)
	pmf <- dbinom(0:N, N,mu/N)
	
	musig <- TBSmeansd(n, pmf, f)
	musig[1] + c(-1,1) * z2 * musig[2]
}

#Exceedance probability for given interval "crit" under alternative:
apower <- function(n, dist, f, crit){
	N <- length(dist)-1
	
	musig <- TBSmeansd(n, dist, f)

	pnorm(crit[1], musig[1], musig[2]) + 1-pnorm(crit[2], musig[1], musig[2])
}

#Asymptotic p value:
apvalue <- function(TBS, N,mu, n, f){
	pmf <- dbinom(0:N, N,mu/N)
	
	musig <- TBSmeansd(n, pmf, f)

	2*(1-pnorm(abs(TBS-musig[1])/musig[2]))
}
#apvalue(TBS, N,mu, n, f)







#Six time series with
N <- 8

for(i in 1:6){ #schools
	xe <- data[data[,5]==i,2]
	mue <- mean(xe)
	p0e <- mean(xe==0)
	BIDe <- N*var(xe)/mue/(N-mue)

	print(c( length(xe), mue, p0e, BIDe ))
}
#124   1.8629032   0.2258065   1.5739563
#127   1.3070866   0.4173228   2.0395892
#136   2.3455882   0.1911765   2.0050656
#132   2.1515152   0.2272727   1.8877596
#155   1.6516129   0.3677419   2.1363055
#123   1.804878   0.300813   2.330074





#Use individual means for test design:

N <- 8
tabI <- c(1.1, 1.25, 1.5)
no.I <- length(tabI)

tab.u <- seq(0.01,0.99,0.01)
no.u <- length(tab.u)

#Find optimal u numerically:
for(i in 1:6){ #schools
	xe <- data[data[,5]==i,2]
	mu <- mean(xe)
	n <- length(xe)

	tab.p <- 1-(N-tabI)/(N-1)*(1-mu/N)
	tab.om <- 1-mu/N/tab.p

for(l in 1:no.I){
	p <- tab.p[l]
	om <- tab.om[l]
	pmf <- dzib(0:N, N,p, om)

	f.optim <- function(u){
		f <- function(x) u^x
		-apower(n, pmf, f, crit.bin(N,mu, n, f))
	}
	
	u.opt <- optim(c(1), f.optim, method="L-BFGS-B", lower=c(0.0001), upper=c(0.9999), control=list(ndeps=c(1e-4)))
	apow <- -u.opt$value
	u.opt <- u.opt$par
	f.opt <- function(x) u.opt^x
	print(round(c(n,mu,tabI[l], u.opt, apow, sTBS(xe, N, f.opt), crit.bin(N,mu, n, f.opt), apvalue(sTBS(xe, N, f.opt), N,mu, n, f.opt)), 7))
} #for I
} #for school
#124   1.8629032   1.10   0.5739544   0.1770861   0.7375239   0.8846353   1.1255682   0.0000134
#124   1.8629032   1.25   0.4399113   0.6956804   0.6530943   0.8318271   1.1859325   0.0000820
#124   1.8629032   1.50   0.2566094   0.9993593   0.5441659   0.7307637   1.3078269   0.0012488
#127   1.3070866   1.10   0.6596562   0.1610303   0.6842673   0.9167849   1.0893674   0
#127   1.3070866   1.25   0.4881690   0.6431080   0.5614455   0.8643764   1.1472126   0
#127   1.3070866   1.50   0.2725410   0.9973401   0.4465081   0.7780967   1.2459727   0.0000022
#136   2.3455882   1.10   0.5129428   0.2175514   0.5599128   0.8602693   1.1538022   0
#136   2.3455882   1.25   0.3880090   0.7945873   0.4762320   0.7983062   1.2271810   0.0000009
#136   2.3455882   1.50   0.2335041   0.9999736   0.3919558   0.6847277   1.3717648   0.0002830
#132   2.1515152   1.10   0.5336577   0.2005380   0.5822493   0.8697907   1.1426121   0
#132   2.1515152   1.25   0.4064692   0.7589153   0.4847349   0.8125199   1.2094083   0.0000002
#132   2.1515152   1.50   0.2401374   0.9999007   0.3788286   0.7037766   1.3446641   0.0000790
#155   1.6516129   1.10   0.5619867   0.1996564   0.5425507   0.8948185   1.1132639   0
#155   1.6516129   1.25   0.4124694   0.7755821   0.4258119   0.8431375   1.1712930   0
#155   1.6516129   1.50   0.2242685   0.9999249   0.3154445   0.7526457   1.2776223   0.0000002
#123   1.8048780   1.10   0.5825806   0.1735953   0.5812968   0.8877901   1.1219784   0
#123   1.8048780   1.25   0.4462225   0.6850562   0.4999434   0.8354506   1.1815776   0
#123   1.8048780   1.50   0.2597058   0.9991581   0.4186713   0.7359955   1.3008356   0.0000315







####################
#Further Simulations
####################

#For comparison, try out binomial U-test,
#see Savage (1970, p. 188)

Ubin <- function(x, N){
	n <- length(x)
	R <- sum(x)
	v <- 2*n * (N*n/(N*n-1))^2 * (n*N-R-1)/(n*N-R) * (R-1)/R * n*N/(n*N-2) * (N-1)/(n*N-3)
	sqrt((n-1)/v) * (var(x)/mean(x)/(1-mean(x)/N) - N*n/(N*n-1))
}
level <- 0.05
crit <- c(-1,1) * qnorm(1-level/2)


set.seed(123)
reps <- 1e5

N <- 8

no.n <- 6
tabmu <- tabn <- rep(NA, no.n)
for(i in 1:no.n){ #schools
	xe <- data[data[,5]==i,2]
	tabmu[i] <- mean(xe)
	tabn[i] <- length(xe)
}
tabI <- c(1.1, 1.25, 1.5)
no.I <- length(tabI)



#Finite-sample performance of sizes:
reject <- array(0, c(no.I, no.n, 2))
for(k in 1:no.n){
	n <- tabn[k]
	mu <- tabmu[k]
	
for(l in 1:no.I){
	Id <- tabI[l]
	p <- 1-(N-Id)/(N-1)*(1-mu/N)
	om <- 1-mu/N/p

	for(r in 1:reps){
		data0 <- rbinom(n, N,mu/N)
		data1 <- rzib(n, N,p, om)
		
		T0 <- Ubin(data0, N)
		T1 <- Ubin(data1, N)
		
		if(T0<crit[1] || T0>crit[2]) reject[l,k,1] <- reject[l,k,1]+1
		if(T1<crit[1] || T1>crit[2]) reject[l,k,2] <- reject[l,k,2]+1
	} #for reps

	#print(c(mu,Id,n, opt.des[l,k,], reject[l,k,]))
}} #for n,I

results <- array(NA, c(no.I*no.n, 6))
i <- 0
for(k in 1:no.n){
for(l in 1:no.I){
	i <- i+1
	results[i,] <- c(N, tabmu[k], tabI[l], tabn[k], reject[l,k,])
}} #for n,I
results

 # [1,]    8 1.862903 1.10  124 4827 14999
 # [2,]    8 1.862903 1.25  124 4955 53803
 # [3,]    8 1.862903 1.50  124 4842 96466
 # [4,]    8 1.307087 1.10  127 4939 14818
 # [5,]    8 1.307087 1.25  127 4902 54050
 # [6,]    8 1.307087 1.50  127 4948 97014
 # [7,]    8 2.345588 1.10  136 4886 16261
 # [8,]    8 2.345588 1.25  136 5018 56992
 # [9,]    8 2.345588 1.50  136 4770 96982
# [10,]    8 2.151515 1.10  132 4920 15849
# [11,]    8 2.151515 1.25  132 4974 56024
# [12,]    8 2.151515 1.50  132 4934 97072
# [13,]    8 1.651613 1.10  155 4939 17055
# [14,]    8 1.651613 1.25  155 4949 61888
# [15,]    8 1.651613 1.50  155 4922 98803
# [16,]    8 1.804878 1.10  123 5026 14742
# [17,]    8 1.804878 1.25  123 5078 53377
# [18,]    8 1.804878 1.50  123 4966 96327








#For comparison, try out zero index of Kim et al. (2018):

ZIbin <- function(x, N){
	n <- length(x)
	m <- mean(x)
	p0 <- mean(x==0)
	
	sqrt(n)*(p0/(1-m/N)^N-1 + (N-1)/2/n * m/N/(1-m/N)) / sqrt( 1/(1-m/N)^N-1 - m/(1-m/N) )
}
level <- 0.05
crit <- c(-1,1) * qnorm(1-level/2)


set.seed(123)
reps <- 1e5

N <- 8

no.n <- 6
tabmu <- tabn <- rep(NA, no.n)
for(i in 1:no.n){ #schools
	xe <- data[data[,5]==i,2]
	tabmu[i] <- mean(xe)
	tabn[i] <- length(xe)
}
tabI <- c(1.1, 1.25, 1.5)
no.I <- length(tabI)



#Finite-sample performance of sizes:
reject <- array(0, c(no.I, no.n, 2))
for(k in 1:no.n){
	n <- tabn[k]
	mu <- tabmu[k]
	
for(l in 1:no.I){
	Id <- tabI[l]
	p <- 1-(N-Id)/(N-1)*(1-mu/N)
	om <- 1-mu/N/p

	for(r in 1:reps){
		data0 <- rbinom(n, N,mu/N)
		data1 <- rzib(n, N,p, om)
		
		T0 <- ZIbin(data0, N)
		T1 <- ZIbin(data1, N)
		
		if(T0<crit[1] || T0>crit[2]) reject[l,k,1] <- reject[l,k,1]+1
		if(T1<crit[1] || T1>crit[2]) reject[l,k,2] <- reject[l,k,2]+1
	} #for reps

	#print(c(mu,Id,n, opt.des[l,k,], reject[l,k,]))
}} #for n,I

results <- array(NA, c(no.I*no.n, 6))
i <- 0
for(k in 1:no.n){
for(l in 1:no.I){
	i <- i+1
	results[i,] <- c(N, tabmu[k], tabI[l], tabn[k], reject[l,k,])
}} #for n,I
results

 # [1,]    8 1.862903 1.10  124 4962 21658
 # [2,]    8 1.862903 1.25  124 5000 75887
 # [3,]    8 1.862903 1.50  124 4936 99770
 # [4,]    8 1.307087 1.10  127 4982 17714
 # [5,]    8 1.307087 1.25  127 5009 68436
 # [6,]    8 1.307087 1.50  127 4964 99583
 # [7,]    8 2.345588 1.10  136 4777 28533
 # [8,]    8 2.345588 1.25  136 4879 85258
 # [9,]    8 2.345588 1.50  136 4735 99948
# [10,]    8 2.151515 1.10  132 4752 25573
# [11,]    8 2.151515 1.25  132 4950 82003
# [12,]    8 2.151515 1.50  132 5014 99910
# [13,]    8 1.651613 1.10  155 5142 23451
# [14,]    8 1.651613 1.25  155 5180 81549
# [15,]    8 1.651613 1.50  155 5016 99965
# [16,]    8 1.804878 1.10  123 5056 20925
# [17,]    8 1.804878 1.25  123 5065 75187
# [18,]    8 1.804878 1.50  123 4951 99779
