# Function that extracts information when running the simulation study
ss.extract <- function(out, ciMat, muMat, sig){
	output <- NULL
	mu.int <-  aperm(apply(out$mu,c(1,2),emp.hpd), c(3,1,2))	
	sig2.int <-  aperm(apply(out$sig2,c(1,2),emp.hpd), c(3,1,2))	

	rho <- list()
	cov.mu <- cov.sig2 <- numeric()
	rand <- rand.cont <- rand.lag <- numeric()
	rand.mat <- matrix(NA, Tm, Tm)
	for(k in 1:Tm){
		epam <- expectedPairwiseAllocationMatrix(t(out$Si[k,,]))
		rho[[k]] <- salso(epam, structure="clustering", loss="binder")
		rand[k] <- adj.rand.index(rho[[k]], ciMat[k,])

		cov.mu[k] <- sum(apply(mu.int[,,k] - muMat[k,],1,prod)	< 0)
		cov.sig2[k] <- sum(apply(sig2.int[,,k] - sig^2,1,prod)	< 0)
	}

	rand.first.last <- adj.rand.index(rho[[1]], rho[[Tm]])

	for(k in 1:(Tm-1)){
		rand.cont[k] <- adj.rand.index(rho[[k]], rho[[k+1]])
		rand.lag[k] <- adj.rand.index(rho[[1]], rho[[k+1]])
	}	

	for(k in 1:Tm){
		for(kk in 1:Tm){
			rand.mat[k,kk] <- adj.rand.index(rho[[k]], rho[[kk]])
		}
	}

	mean.cov.mu <- sum(cov.mu)/(N*Tm)
	mean.cov.sig2 <- sum(cov.sig2)/(N*Tm)
	output$rand <- rand
	output$rand.cont <- rand.cont
	output$rand.lag <- rand.lag
	output$rand.first.last <- rand.first.last
	output$mean.cov.mu <- mean.cov.mu
	output$mean.cov.sig2 <- mean.cov.sig2
	output$rand.mat <- rand.mat
	output
}



# Function that computes ``robust'' lpmls and waic
lpml.robust <- function(llike){

	omega <- 1/exp(llike)

	lpml <- sum(log(1/(apply(omega,2,mean))))

	omegabar <- apply(omega,2,mean)

	omegatil <- matrix(NA, nrow=nrow(omega), ncol=ncol(omega))
	for(i in 1:nrow(omega)){
		for(j in 1:ncol(omega)){

			omegatil[i,j] <- min(omega[i,j], sqrt(nrow(omega)*omegabar[j]))
		}
	}

	lpml.robust <- sum(log(1/(apply(omegatil,2,mean))))

	# Compute lpml by droping iteration that produces inf
	indx <- which(llike < -709, arr.ind=TRUE)[,1]
	CPO <- apply(llike[-indx,], 2, function(x) mean(1/exp(x)))
	if(length(indx) == 0){
		CPO <- apply(llike, 2, function(x) mean(1/exp(x)))
	}
	lpml.drop.inf <- sum(log(1/CPO))

	# Compute lpml by removing all llike values less than
    # the 5th quantile
    mn995 <- function(x) mean(x[x < quantile(x, 0.995)])

	lpml.quant <- sum(log(1/(apply(omega,2,mn995))))


	# Comput lpml by summing off min (see paper for details)
	nobs <- ncol(llike)
	niter <- nrow(llike)
	minll <- apply(llike,2,min)
	tmp <- numeric()
	for(i in 1:nobs){

		tmp[i] <- sum(exp(minll[i] - llike[,i]))

	}

	lpml.log.scale <- sum(log(niter) + minll - log(tmp))

	# compute the waic
  mnllike <- apply(llike, 2, mean)
  mnlike <- apply(exp(llike), 2, mean)
  elppWAIC <- sum(2*mnllike - log(mnlike))
  waic <- -2*elppWAIC

	out <- c(lpml, lpml.robust, lpml.log.scale, lpml.drop.inf, lpml.quant, waic)
	names(out) <- c("lpml",
                    "lpml.stable.gelman",
                    "lpml.stable.log",
                    "lpml.drop.inf",
                    "lpml.drop.quant995",
                    "waic")
	out
}






## Marginal likelihood for auxilliary and double dipper cohesion functions
#sk <- s_std[1:75,]; mu0<- apply(s_std,2,mean);k0<-1;v0<-2;L0<- diag(2)
marglik <- function(sk, mu0, k0, v0, L0,log=TRUE){
	G2a <- function(a,log=TRUE) {
		out <- 0.5*log(pi) + lgamma(a) + lgamma(a-0.5)
		if(!log) out <- exp(out)
		out
	}

#	cat("sk = ", sk, "\n")

	n <- dim(sk)[1]
	sbar <- apply(sk,2,mean)
    Vs <- matrix(0,2,2)
    for(jj in 1:n){
    	Vs <- Vs + (sk[jj,] - sbar) %*% t(sk[jj,] - sbar)
    }
	mun <- k0/(k0+n)*mu0 + n/(n+k0)*sbar
	kn <- k0+n
	vn <- v0+n
	Ln = L0 + Vs + k0*n/(k0+n)*(sbar-mu0)%*%t(sbar-mu0)
	
	munn <- kn/(kn+n)*mun + n/(n+kn)*sbar
	knn <- kn + n
	vnn <- vn + n
	Lnn <- Ln + Vs + kn*n/(kn+n)*(sbar-mun)%*%t(sbar-mun)

#	cat("det(L0)^(0.5*v0)/det(Ln)^(0.5*vn) = ", det(L0)^(0.5*v0)/det(Ln)^(0.5*vn), "\n")

	aux <- -n*log(pi) + (G2a(0.5*vn,log=TRUE) - G2a(0.5*v0,log=TRUE)) + 
	                     ((0.5*v0)*log(det(L0)) - (0.5*vn)*log(det(Ln))) + 
	                     (0.5*2)*log(k0/kn)
	
	dd <- (-n)*log(pi) + (G2a(0.5*vnn,log=TRUE) - G2a(0.5*vn,log=TRUE)) + 
	                     ((0.5*vn)*log(det(Ln)) - (0.5*vnn)*log(det(Lnn))) + 
	                     (0.5*2)*log(kn/knn)

	if(!log){aux<-exp(aux); dd<-exp(dd)}

	c(aux,dd)

}	

	



## I want to draw a cluster configuration from the prior
## Use the prediction probabilies (equation (8) in Quintana)?
## n <- 1; cohesion <- 4; M <- 1; aq <- 0.5; a1 <-1; plot<-FALSE
rPPMs <- function(n, s, cohesion, M, aq,a1=1, parm=c(0,1,2,1), ppart=NULL, plot=FALSE){
	
	out <- NULL
	S <- NULL
	ndraw <- n

	colnames(s) <- c("s1", "s2")#, "x")
	sbar <- apply(s,2,mean)
	ssd <- apply(s,2,sd)
	s_std <- t((t(s) - sbar)/ssd)

	a <- quantile(dist(s), probs=aq)

	m <- rep(C3_4parm[1], 2)
	k0 <- C3_4parm[2]
	nu0 <- C3_4parm[3]
	L0 <- C3_4parm[4]*diag(2)

	# I only want to assign cluster labels for subjects 
	# that are not a part of ppart.
	
	nloc <- dim(s)[1]
	
	
##	upord <- 1:nloc
	## Try to create a simulation study to investigate properties of prior a bit more.
	for(j in 1:ndraw){
		
		upord <- (sample(1:nloc, nloc, replace=FALSE))
		if(!is.null(ppart))	upord <- sample(c(1:nloc)[ppart==0], sum(ppart==0), replace=FALSE)

		if(plot){
			plot(s, type='n')
			points(s[upord[1],,drop=FALSE],pch="a", cex=0.85)
		}	
		if(is.null(ppart)){
			Si <- rep(0,nloc)
			Si[upord[1]] <- 1
			upord <- upord[-1]
			## vector indicating number of obs in each cluster
			nh <- 1 

			## this is the number of clusters
			nk <- 1


		} else {
			Si <- ppart

			## vector indicating number of obs in each cluster
			nh <- tabulate(Si)

			## this is the number of clusters
			nk <- length(nh)

		}
		
		## vector to story probabilities associated with each cluster assigned to evaluate density
		probkeep <- rep(1,nloc)
		for(i in upord){
		
			ph <- numeric()

		
			for(k in 1:nk){
				stmp0 <- rbind(s_std[Si==k,,drop=FALSE]) 
				stmp <- rbind(s_std[Si==k,,drop=FALSE],s_std[i,])
	
#				print(stmp0)
#				print(stmp)
				
	
				pwd0 <- dist(stmp0)
				pwd1 <- dist(stmp)


				smn0 <- matrix(apply(stmp0,2, mean), nrow=1, byrow=TRUE)
				dif0 <- (stmp0 - smn0[rep(1, each=nh[k]),])^2
				ssq0 <- sqrt(apply(dif0,1,sum))
				dist0 <- sum(ssq0) 		

				smn1 <- matrix(apply(stmp,2, mean), nrow=1, byrow=TRUE)
				dif1 <- (stmp - smn1[rep(1, each=nh[k]+1),])^2
				ssq1 <- sqrt(apply(dif1,1,sum)) ## These are all distances to centroid
				dist1 <- sum(ssq1) ## These are the sum of all centroid distances		

				## Recall that all cohesions treat singletons differently and that is why
				## the denominator depends of this ratio depends on cluster size (by construction
				## the numerator cluster size is always greater than 1).

				## Care must be taken with cluster singletons.  This only happens here in the denominator or if
				## considering assiging to new cluster.

				## I change the cohesion 1 to explicitly depend on whether cluster is singleton or not.  If it is a 
				## singleton the cohesion simple returns M.  The remaining cohesions internally take care of the case
				## of singletons
				
				if(cohesion==1){			

					lnum <- log(M) + lgamma(nh[k]+1)-
					        (lgamma(a1*dist1)*(dist1>=1) + log(dist1)*(dist1 < 1))
					if(nh[k] > 1){
						lden <- log(M) + lgamma(nh[k])-
						        (lgamma(a1*dist0)*(dist0>=1) + log(dist0)*(dist0<1))
					}else{
						lden <- log(M)
					}
					lC0 <- log(M)
				}
				
				if(cohesion==2){
#					lnum <- log(M)+lgamma(nh[k]+1) - lfactorial(1+sum(pwd1))
#					if(nh[k] >1)lden<-log(M)+lgamma(nh[k])-lfactorial(1+sum(pwd0))else lden<-log(M)
					lnum <- log(M) + lgamma(nh[k]+1) + log(prod(pwd1<a)==1)
					if(nh[k] > 1){
						lden <- log(M) + lgamma(nh[k]) + log(prod(pwd0<a)==1)
					}else{
						lden <- log(M)	
					} 
					lC0 <- log(M)
				}


				if(cohesion==3){
					Llik0 <- marglik(stmp0, m, k0, nu0, L0,log=TRUE)[1]
					Llik <- marglik(stmp, m, k0, nu0, L0,log=TRUE)[1]

					lnum <- log(M)+lgamma(nh[k]+1) + Llik
					lden <- log(M)+lgamma(nh[k]) + Llik0 
					lC0 <- log(M) + marglik(s_std[i,,drop=FALSE], m, k0, nu0, L0,log=TRUE)[1]
				}

				if(cohesion==4){
#					print(stmp)
#					print(stmp0)
					Llik0 <- marglik(stmp0, m, k0, nu0, L0,log=TRUE)[2]
					Llik <- marglik(stmp, m, k0, nu0, L0,log=TRUE)[2]

#					cat("Llik0 = ", Llik0, "\n")
#					cat("Llik = ", Llik, "\n")

					lnum <- log(M)+lgamma(nh[k]+1) + Llik
					lden <- log(M)+lgamma(nh[k]) + Llik0 

#					cat("lnum = ", lnum, "\n")
#					cat("lden = ", lden, "\n")

#					lnum <- log(M)+lgamma(nh[k]+1) 
#					lden <- log(M)+lgamma(nh[k])  
#					print(s_std[i,,drop=FALSE])
					lC0 <- log(M)+ marglik(s_std[i,,drop=FALSE], m, k0, nu0, L0,log=TRUE)[2]
#					cat("lC0 = ", lC0, "\n")
#					lC0 <- log(M)
				}
				
				
#				if(!cohesion%in%c(2,3)) ph[k] <- num/den
#				if(cohesion%in%c(2,3)) ph[k] <- lnum - lden
#				if(!cohesion%in%c(3,4)) ph[k] <- num/den
#				if(cohesion%in%c(3,4)) ph[k] <- lnum - lden
#				cat("lnum = ", lnum, "\n")
#				cat("lden = ", lden, "\n")
				ph[k] <- lnum - lden		
			}
#			cat("k = ", k, "\n")
#			cat("nk = ", nk, "\n")
			ph[nk+1] <- lC0
#			cat("ph = ", round(ph,3), "\n")

			ph <- exp(ph - max(ph))
			
			probs <- ph/sum(ph)

#			cat("probs = ", round(probs,3), "\n")
			tmp <- sample(1:(nk+1), size=1, prob=probs)
			Si[i] <- tmp
			nh <- tabulate(Si)
			nk <- length(nh)
#			print(Si)
#			print(nh)
#			print(nk)
			probkeep[i] <- probs[tmp] ## This is to density associated with 
			if(plot){
				Sys.sleep(0.05)
				cat("point", i, "probs = ", probs, "\n")
				points(s[i,,drop=FALSE],pch=letters[Si[i]], cex=0.85, col=Si[i])
			}
		}
		Si <- as.numeric(as.character(factor(Si, labels=order(unique(Si), decreasing=FALSE))))

		S <- rbind(S,Si)

	}

	out$S <- S
	out$ldensity <- sum(log(probkeep))
	out

}






## Functions needed to generate data

# November 2 2018
# N - number of experimental units
# M - DP scale (dispesrsion) parameter
# rho - temporal dependence parameter
# ntime - number of time points
# tau - sd associated with atom generation
# sig - sd associated with data generation
# Caron - logical indicating if temporal dependence follows Caron
# FirstPart - Allows me to provide the first partition if so desired
# Space - logical indicating of sPPM should be used to generate partition or not
# N <- 3;M <- 1;rho <- 0.5;ntime <- 2;Caron=FALSE;FirstPart=NULL
rtpartition <- function(N,M,rho,ntime,Caron=TRUE, FirstPart=NULL){
	library(plyr)
	out <- NULL

	ci <- FirstPart

	if(is.null(FirstPart)){
		# Create cluster at time period t = 1
		ci <- 1 # cluster label
		mh <- 1 # vector of cluster sizes
		K <- 1 # number of clusters
		Tm <- ntime

		for(k in 2:N){

			p <- c(mh/(k-1+M), M/(k-1+M))
			ci <- c(ci, sample(1:(K+1), 1, prob=p))
			mh <- table(ci)
			K <- length(unique(ci))
		}

	}

#	print(ci)
#	print(mustar)
#	plot(1:N, ci, col=ci, pch=19, cex=1.25)

	# Now to generate 10-1=9 more partitions that are temporally dependent
	# Will first try uniform deletion
	ciMat <- matrix(NA, nrow=Tm, ncol=N)
	ciMat[1,]  <- ci
	K <-  length(unique(ci))
	gMat <- matrix(0,nrow=Tm, ncol=N)

#	print(ciMat)
	if(Tm > 1){
		if(Caron){
			for(t in 2:Tm){
				cit <- NULL
				r <- rbinom(1, sum(mh), 1-rho)
				dnk <- sample(1:length(ci), r)
				if(r > 0){
					ci <- ci[-dnk]
				}	
				mh <- tabulate(ci)
#				K <- length(unique(ci))
				K <- length(mh)

				for(k in 1:N){

					p <- c(mh/(sum(mh)+M), M/(sum(mh)+M))
					cit[k] <- sample(1:(K+1), 1, prob=p)
					ci <- c(ci, cit[k])
					mh <- table(ci)
					K <- length(unique(ci))
				}
	
				ciMat[t,] <- cit
	
			}

		}else{
		
			for(t in 2:Tm){
#				cat("t = ", t, "\n")
				r <- rbinom(1, sum(mh), 1-rho)
				dnk <- sample(1:length(ci), r)
				if(r > 0){
					ci[dnk]<- 0
					gMat[t,dnk] <- 1

				}	

#				print(ci)

				mh <- tabulate(ci[ci!=0]);
				if((K - length(mh)) > 0){
					mh <- c(mh, rep(0, K-length(mh)))
				}

				K <- length(unique(ci[ci!=0]))


#				print(K); print(length(mh))

#				print(mh)

				if(r < N){
				
					ci[ci!=0]<-as.numeric(as.character(factor(ci[ci!=0], labels=order(unique(ci[ci!=0]),decreasing=FALSE))))

				}
#				print(ci)



				for(k in dnk){
#					cat("k = ", k, "\n")
#					print(mh)
					p <- 1
					if(K>0) p <- c(mh[mh!=0]/(sum(mh[mh!=0])+M), M/(sum(mh[mh!=0])+M))
#					print(p)
#					cat("K + 1 = ", K+1, "\n")

					ci[k] <- sample(1:(K+1), 1, prob=p)
					mh <- table(ci[ci!=0])
					K <- length(unique(ci[ci!=0]))
#					print(ci)
#					print(table(ci))
				}
#				print(ci)
				citmp <- factor(ci,  labels=order(unique(ci), decreasing=FALSE))

#				print(citmp)
				ciMat[t,] <- as.numeric(as.character(citmp))
#				print(ciMat)

			}
		
		}
	}
	out$ciMat <- ciMat
	out$gMat <- gMat
	out
}


## Functions needed to generate data

# November 2 2018
# N - number of experimental units
# M - DP scale (dispesrsion) parameter
# rho - temporal dependence parameter
# ntime - number of time points
# tau - sd associated with atom generation
# sig - sd associated with data generation
# Caron - logical indicating if temporal dependence follows Caron
# FirstPart - Allows me to provide the first partition if so desired
# Space - logical indicating of sPPM should be used to generate partition or not
# N <- 225;M <- 1;rho <- 0.5;ntime <- 5;Caron=FALSE;FirstPart=NULL
# cohesion = 3; C3_4parm=c(0,1,5,1)
rtpartition.space <- function(N,M,rho,ntime,Caron=TRUE, FirstPart=NULL,s_coords=NULL,cohesion=NULL,C3_4parm=NULL){

	out <- NULL

	ci <- FirstPart
	Tm <- ntime

	if(is.null(FirstPart)){
		# Create cluster at time period t = 1
		ci <- 1 # cluster label
		mh <- 1 # vector of cluster sizes
		K <- 1 # number of clusters

		for(k in 2:N){

			p <- c(mh/(k-1+M), M/(k-1+M))
			ci <- c(ci, sample(1:(K+1), 1, prob=p))
			mh <- table(ci)
			K <- length(unique(ci))
		}

		if(!is.null(s_coords)){
			ci <- rPPMs(1, s_coords, cohesion, M, aq=0.95,a1=1, C3_4parm=C3_4parm, ppart=NULL, plot=FALSE)[[1]][1,]
			mh <- table(ci)
			K <- length(unique(ci))
		}
	}

#	print(ci)
#	print(mustar)
#	plot(1:N, ci, col=ci, pch=19, cex=1.25)
	ci <- as.numeric(as.character(factor(ci,  labels=order(unique(ci), decreasing=FALSE))))
#	print(citmp)
	# Now to generate Tm-1 more partitions that are temporally dependent
	# Will first try uniform deletion
	ciMat <- matrix(NA, nrow=Tm, ncol=N)
	ciMat[1,]  <- ci
	K <-  length(unique(ci))
	gMat <- matrix(0,nrow=Tm, ncol=N)

#	print(ciMat)
	if(Tm > 1){

		if(Caron){
			for(t in 2:Tm){
				cit <- NULL
				r <- rbinom(1, sum(mh), 1-rho)
				dnk <- sample(1:length(ci), r)
				if(r > 0){
					ci <- ci[-dnk]
				}	
				mh <- tabulate(ci)
#				K <- length(unique(ci))
				K <- length(mh)

				for(k in 1:N){

					p <- c(mh/(sum(mh)+M), M/(sum(mh)+M))
					cit[k] <- sample(1:(K+1), 1, prob=p)
					ci <- c(ci, cit[k])
					mh <- table(ci)
					K <- length(unique(ci))
				}
	
				ciMat[t,] <- cit
	
			}

		}else{
			for(t in 2:Tm){
#				cat("t = ", t, "\n")
				r <- rbinom(1, sum(mh), 1-rho)
				dnk <- sample(1:length(ci), r)
				if(r > 0){
					ci[dnk]<- 0
					gMat[t,dnk] <- 1
	
				}	

#				print(ci)

				mh <- tabulate(ci[ci!=0]);
				if((K - length(mh)) > 0){
					mh <- c(mh, rep(0, K-length(mh)))
				}

				K <- length(unique(ci[ci!=0]))


#				print(K); print(length(mh))

#				print(mh)
				# need to relabel cluster labesl associated with those that
				# are not being reallocated.
				if(r < N){
					ci[ci!=0] <- as.numeric(as.character(factor(ci[ci!=0],  labels=order(unique(ci[ci!=0]), decreasing=FALSE))))
				}
#				print(ci)



				if(is.null(s_coords)){	
					for(k in dnk){
#						cat("k = ", k, "\n")
#						print(mh)
						p <- 1
						if(K>0) p <- c(mh[mh!=0]/(sum(mh[mh!=0])+M), M/(sum(mh[mh!=0])+M))
#						print(p)
#						cat("K + 1 = ", K+1, "\n")

						ci[k] <- sample(1:(K+1), 1, prob=p)
						mh <- table(ci[ci!=0])
						K <- length(unique(ci[ci!=0]))
#						print(ci)
#						print(table(ci))
					}
			
				}
#				print(ci)
				if(!is.null(s_coords)){
					ppart <- ci
					if(sum(ci)==0) ppart <- NULL
					ci <- rPPMs(1, s_coords, cohesion, M, aq=0.95,a1=1, C3_4parm=C3_4parm, ppart=ppart, plot=FALSE)[[1]][1,]

				}

#				print(ci)
				citmp <- factor(ci,  labels=order(unique(ci), decreasing=FALSE))

#				print(citmp)



				ciMat[t,] <- as.numeric(as.character(citmp))
			}
		
		}
	}
#	print(ciMat)
	out$ciMat <- ciMat
	out$gMat <- gMat
	out
}

# phi1 = 0, mus generated randomly
datgenFULL <- function(N,M,rho,ntime,Caron=TRUE, FirstPart=NULL, tau=5,sig=1,eta1=0,
							phi0=0,phi1=0, lam=1){
	library(rmutil, warn.conflicts=FALSE)
	out <- NULL
	Tm <- ntime
	muMat <- matrix(0, nrow=Tm, ncol=N)
	YMat <- matrix(0, nrow=Tm, ncol=N)

	# Generate partition
	out <- rtpartition(N=N, M=M, rho=rho, ntime=ntime,Caron=Caron, FirstPart=FirstPart)
	ciMat <- out$ciMat
	gMat <- out$gMat

	# Generate thetas
	theta <- numeric()
	theta[1] <- rnorm(1, phi0, lam)
	for(t in 2:ntime){
		theta[t] <- rnorm(1, phi0*(1-phi1) + phi1*theta[t-1], lam*sqrt(1-phi1^2))
	}

	# Generate mustar and sigstar
	for(t in 1:ntime){
		K <- length(unique(ciMat[t,]))
		mustar <- rnorm(K, theta[t], tau)
		muMat[t,] <- mustar[ciMat[t,]]

	}

	# Generate xi and eta1i if eta1 is not specified
	if(is.null(eta1)){
		xi <- rlaplace(N, 0, 1)
		eta1 <- (exp(xi) - 1)/(1 + exp(xi))
	}
	
	# Generate observations
	YMat[1,] <- rnorm(N, muMat[1,], sig)
	for(t in 2:ntime){
		
		YMat[t,] <- rnorm(N, muMat[t,] + eta1*YMat[t-1,], sig*sqrt(1-eta1^2))
	}
	
	out$ciMat <- ciMat
	out$gMat <- gMat
	out$theta <- theta
	out$eta1 <- eta1
	out$muMat <- muMat
	out$YMat <- YMat
	out	

}
# tmp <- datgenFULL(N=50, M=1, rho=0.75, ntime=5, Caron=FALSE, FirstPart=NULL, 
#					tau=5, sig=1, eta1=0.75, phi0=0, phi1=0.75, lam=1)



# phi1 = 0, mus generated randomly
datgenFULLspace <- function(N,M,rho,ntime,Caron=TRUE, FirstPart=NULL, tau=5,sig=1,eta1=0,
							phi0=0,phi1=0,lam=1,s_coords=NULL,cohesion=NULL, C3_4parm=NULL){

	library(rmutil, warn.conflicts=FALSE)
	out <- NULL
	Tm <- ntime
	muMat <- matrix(0, nrow=Tm, ncol=N)
	YMat <- matrix(0, nrow=Tm, ncol=N)

	# Generate partition
	out <- rtpartition.space(N=N, M=M, rho=rho, ntime=ntime,Caron=Caron,FirstPart=FirstPart,
							s_coords=s_coords,cohesion=cohesion, C3_4parm=C3_4parm )
	ciMat <- out$ciMat
	gMat <- out$gMat

	# Generate thetas
	theta <- numeric()
	theta[1] <- rnorm(1, phi0, lam)
	for(t in 2:ntime){
		theta[t] <- rnorm(1, phi0*(1-phi1) + phi1*theta[t-1], lam*sqrt(1-phi1^2))
	}

	# Generate mustar and sigstar
	for(t in 1:ntime){
		K <- length(unique(ciMat[t,]))
		mustar <- rnorm(K, theta[t], tau)
		muMat[t,] <- mustar[ciMat[t,]]

	}

	# Generate xi and eta1i if eta1 is not specified
	if(is.null(eta1)){
		xi <- rlaplace(N, 0, 1)
		eta1 <- (exp(xi) - 1)/(1 + exp(xi))
	}
	
	# Generate observations
	YMat[1,] <- rnorm(N, muMat[1,], sig)
	for(t in 2:ntime){
		
		YMat[t,] <- rnorm(N, muMat[t,] + eta1*YMat[t-1,], sig*sqrt(1-eta1^2))
	}
	
	out$ciMat <- ciMat
	out$gMat <- gMat
	out$theta <- theta
	out$eta1 <- eta1
	out$muMat <- muMat
	out$YMat <- YMat
	out	

}


space.time.dat.gen <- function(Tm=2, s_coords){
	
	library("CompRandFld")

	out <- NULL
	# Define the spatial-coordinates of the points:

	N <- nrow(s_coords)
	cmat <- Covmatrix(coordx=s_coords, coordt=c(1:Tm), corrmodel="exp_exp", 
						param=list(scale_s=0.3,scale_t=2,sill=1.75,mean=0))$covmatrix
	stdnorm <- rnorm(nrow(cmat))
	ccmat <- chol(cmat)
	Y <- t(ccmat) %*% stdnorm 
	YMat <- matrix(Y, nrow=N, byrow=TRUE)
	out$YMat <- YMat
	out
	

}

#geodat1 <- list(coords=s_coords, data=YMat[,1])
#geodat2 <- list(coords=s_coords, data=YMat[,2])
#geodat3 <- list(coords=s_coords, data=YMat[,3])

#plot(variog(geodat1))
#variofit(variog(geodat1),cov.model="exponential")
#lines(variofit(variog(geodat1),cov.model="exponential"), lty=3, col='blue')

#par(mfrow=c(2,2))
#image.plot(matrix(YMat[,1], 15, byrow=TRUE))
#image.plot(matrix(YMat[,2], 15, byrow=TRUE))
#image.plot(matrix(YMat[,3], 15, byrow=TRUE))
#image.plot(matrix(YMat[,4], 15, byrow=TRUE))




# This function was used to generate data when exploring correlation on the 
# data level.
# March 17 2018
# N - number of experimental units
# M - DP scale (dispesrsion) parameter
# rho - temporal dependence parameter
# ntime - number of time points
# tau - sd associated with atom generation
# sig - sd associated with data generation
# Caron - logical indicating if temporal dependence follows Caron
# FirstPart - Allows me to provide the first partition if so desired
# Type - type of model for the cluster-specific paramters
#		Random - cluster-specific parameters are drawn randomly across time
#		Deterministic - cluster-specific parameters are same over time
#		AR1 - cluster-specific parameters are drawn from an AR(1) process
#				If this is selected, phi0 and phi1 must be supplied.
# N <- 50;M <- 1;rho <- 0.5;ntime <- 5;Caron=FALSE;FirstPart=NULL; tau<-5;sig <- 1;phi0<-0;phi1<-1
rtpartition1 <- function(N,M,rho,ntime,tau,sig,Caron=TRUE, FirstPart=NULL, TYPE="Random",phi0=0,phi1=0){
	out <- NULL
	ci <- FirstPart
	if(is.null(FirstPart)){
		# Create cluster at time period t = 1
		ci <- 1 # cluster label
		mh <- 1 # vector of cluster sizes
		K <- 1 # number of clusters
		Tm <- ntime

		for(k in 2:N){

			p <- c(mh/(k-1+M), M/(k-1+M))
			ci <- c(ci, sample(1:(K+1), 1, prob=p))
			mh <- table(ci)
			K <- length(unique(ci))
		}
	}
	mustar <- rnorm(K, 0, tau)
	if(TYPE == "AR1") mustar <- rnorm(K, phi0, tau)
	mus <- mustar[ci]
#	print(ci)
#	print(mustar)
#	plot(1:N, ci, col=ci, pch=19, cex=1.25)
	# Now to generate 10-1=9 more partitions that are temporally dependent
	# Will first try uniform deletion
	ciMat <- matrix(NA, nrow=Tm, ncol=N)
	ciMat[1,]  <- ci
	K <-  length(unique(ci))
	gMat <- matrix(0,nrow=Tm, ncol=N)
	muMat <- matrix(0, nrow=Tm, ncol=N)
	muMat[1,] <- mus
	YMat <- matrix(NA, nrow=Tm, ncol=N)
	YMat[1,] <- rnorm(N, muMat[1,], sig)
	if(Caron){
		for(t in 2:Tm){
			cit <- NULL
			r <- rbinom(1, sum(mh), 1-rho)
			dnk <- sample(1:length(ci), r)
			if(r > 0){
				ci <- ci[-dnk]
			}	
			mh <- tabulate(ci)
#			K <- length(unique(ci))
			K <- length(mh)

			for(k in 1:N){

				p <- c(mh/(sum(mh)+M), M/(sum(mh)+M))
				cit[k] <- sample(1:(K+1), 1, prob=p)
				ci <- c(ci, cit[k])
				mh <- table(ci)
				K <- length(unique(ci))
			}
	
			ciMat[t,] <- cit
	
		}
	}else{
		
		for(t in 2:Tm){
#			cat("t = ", t, "\n")
			r <- rbinom(1, sum(mh), 1-rho)
			dnk <- sample(1:length(ci), r)
			if(r > 0){
				ci[dnk]<- 0
				gMat[t,dnk] <- 1

			}	

#			print(ci)

			mh <- tabulate(ci[ci!=0]);
			if((K - length(mh)) > 0){
				mh <- c(mh, rep(0, K-length(mh)))
			}

			K <- length(unique(ci[ci!=0]))


#			print(K); print(length(mh))

#			print(mh)
			mustar <- mustar[mh!=0]
#			print(mustar)

			if(r < N){
				
				ci[ci!=0] <- as.numeric(factor(ci[ci!=0], labels=1:length(unique(ci[ci!=0]))))
			}
#			print(ci)



				
			for(k in dnk){
#				cat("k = ", k, "\n")
#				print(mh)
				p <- 1
				if(K>0) p <- c(mh[mh!=0]/(sum(mh[mh!=0])+M), M/(sum(mh[mh!=0])+M))
#				print(p)
#				cat("K + 1 = ", K+1, "\n")
				ci[k] <- sample(1:(K+1), 1, prob=p)
				mh <- table(ci[ci!=0])
				K <- length(unique(ci[ci!=0]))
#				print(ci)
#				print(table(ci))
			}
			if(TYPE=="Random") mustar <- rnorm(length(mustar),0,tau)
			if(TYPE=="Deterministic") mustar <- mustar
			if(TYPE=="AR1") mustar <- rnorm(length(mustar), phi0+mustar*phi1, tau*sqrt(1-phi1^2))
			
#			cat("K=", K, "\n")
#			print(mustar)	
			# Need to add mustar for any new clusters
			mustar <- c(mustar, rnorm(K - length(mustar), 0, tau))
			if(TYPE == "AR1") mustar <- c(mustar, rnorm(K - length(mustar), phi0, tau))
#			print(mustar)

			ciMat[t,] <- as.numeric(factor(ci, labels=1:length(unique(ci))))
			muMat[t,] <- mustar[ci]

			YMat[t,] <- rnorm(N, muMat[t,], sig)
			
		}
			
		
	}
	
	out$ciMat <- ciMat
	out$gMat <- gMat
	out$muMat <- muMat
	out$YMat <- YMat
	out
}







mattapply <- function(mat, INDEX, FUN){
	
	dms <- dim(mat)
	out <- NULL
	for(i in 1:dms[2]){
				
		out <- cbind(out,tapply(mat[,i], INDEX, FUN))
				
	}
	
	out			

}
	


## Function that creates a color scale used in spatial plots
myColorRamp <- function(colors, values) {
    v <- (values - min(values))/diff(range(values))
    x <- colorRamp(colors)(v)
    rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}
