# This file explores procedure outline in 
# Generalized Polya Urn for Time-varying Dirichlet Process Mixtures
# by Caron et al (see literature review folder)
# 
# We generate a sequence of partitions and see how they
# evolve over discrete time based.

# Set working directory to location of JCGS_Codes
setwd("JCGS_codes")
source("Functions.R")

library(MixSim)
library(mclust)
library(fields)
out.ddp <- list()

out.caron005 <- list()
out.caron025 <- list()
out.caron050 <- list()
out.caron075 <- list()
out.caron095 <- list()
out.caron100 <- list()

out.ours005 <- list()
out.ours025 <- list()
out.ours050 <- list()
out.ours075 <- list()
out.ours095 <- list()
out.ours100 <- list()
# Generate partitions to see differences between our approach and that of Caron
N <- 20
M <- 0.5
nsim <- 10000
nt <- 10
omat1 <- omat2 <- omat3 <- omat4 <- omat5 <- omat6 <- omat7 <- omat8 <- omat9 <- omat10 <- omat11 <- omat12 <- matrix(NA, nt, nt)
set.seed(1)
for(ii in 1:nsim){
	if(ii %% 1000 == 0) cat("ii = ", ii, "\n")
	caron005 <- rtpartition1(N=N,M=M,rho=0.05,ntime=nt,tau=1,sig=0.1,Caron=TRUE, FirstPart=NULL, TYPE="Random",phi0=0,phi1=0)
	caron025 <- rtpartition1(N=N,M=M,rho=0.25,ntime=nt,tau=1,sig=0.1,Caron=TRUE, FirstPart=NULL, TYPE="Random",phi0=0,phi1=0)
	caron050 <- rtpartition1(N=N,M=M,rho=0.50,ntime=nt,tau=1,sig=0.1,Caron=TRUE, FirstPart=NULL, TYPE="Random",phi0=0,phi1=0)
	caron075 <- rtpartition1(N=N,M=M,rho=0.75,ntime=nt,tau=1,sig=0.1,Caron=TRUE, FirstPart=NULL, TYPE="Random",phi0=0,phi1=0)
	caron095 <- rtpartition1(N=N,M=M,rho=0.95,ntime=nt,tau=1,sig=0.1,Caron=TRUE, FirstPart=NULL, TYPE="Random",phi0=0,phi1=0)
	caron100 <- rtpartition1(N=N,M=M,rho=1,ntime=nt,tau=1,sig=0.1,Caron=TRUE, FirstPart=NULL, TYPE="Random",phi0=0,phi1=0)

	ours005 <- rtpartition1(N=N,M=M,rho=0.05,ntime=nt,tau=1,sig=0.1,Caron=FALSE, FirstPart=NULL, TYPE="Random",phi0=0,phi1=0)
	ours025 <- rtpartition1(N=N,M=M,rho=0.25,ntime=nt,tau=1,sig=0.1,Caron=FALSE, FirstPart=NULL, TYPE="Random",phi0=0,phi1=0)
	ours050 <- rtpartition1(N=N,M=M,rho=0.50,ntime=nt,tau=1,sig=0.1,Caron=FALSE, FirstPart=NULL, TYPE="Random",phi0=0,phi1=0)
	ours075 <- rtpartition1(N=N,M=M,rho=0.75,ntime=nt,tau=1,sig=0.1,Caron=FALSE, FirstPart=NULL, TYPE="Random",phi0=0,phi1=0)
	ours095 <- rtpartition1(N=N,M=M,rho=0.95,ntime=nt,tau=1,sig=0.1,Caron=FALSE, FirstPart=NULL, TYPE="Random",phi0=0,phi1=0)
	ours100 <- rtpartition1(N=N,M=M,rho=1,ntime=nt,tau=1,sig=0.1,Caron=FALSE, FirstPart=NULL, TYPE="Random",phi0=0,phi1=0)

	for(j in 1:nt){
		for(jj in 1:nt){
			omat1[j,jj] <- adjustedRandIndex(caron005$ciMat[j,], caron005$ciMat[jj,])			
			omat2[j,jj] <- adjustedRandIndex(caron025$ciMat[j,], caron025$ciMat[jj,])			
			omat3[j,jj] <- adjustedRandIndex(caron050$ciMat[j,], caron050$ciMat[jj,])			
			omat4[j,jj] <- adjustedRandIndex(caron075$ciMat[j,], caron075$ciMat[jj,])			
			omat5[j,jj] <- adjustedRandIndex(caron095$ciMat[j,], caron095$ciMat[jj,])			
			omat6[j,jj] <- adjustedRandIndex(caron100$ciMat[j,], caron100$ciMat[jj,])			

			omat7[j,jj] <- adjustedRandIndex(ours005$ciMat[j,], ours005$ciMat[jj,])			
			omat8[j,jj] <- adjustedRandIndex(ours025$ciMat[j,], ours025$ciMat[jj,])			
			omat9[j,jj] <- adjustedRandIndex(ours050$ciMat[j,], ours050$ciMat[jj,])			
			omat10[j,jj] <- adjustedRandIndex(ours075$ciMat[j,], ours075$ciMat[jj,])			
			omat11[j,jj] <- adjustedRandIndex(ours095$ciMat[j,], ours095$ciMat[jj,])			
			omat12[j,jj] <- adjustedRandIndex(ours100$ciMat[j,], ours100$ciMat[jj,])			
		}
	}

	out.caron005[[ii]] <- omat1
	out.caron025[[ii]] <- omat2
	out.caron050[[ii]] <- omat3
	out.caron075[[ii]] <- omat4
	out.caron095[[ii]] <- omat5
	out.caron100[[ii]] <- omat6

	out.ours005[[ii]] <- omat7
	out.ours025[[ii]] <- omat8
	out.ours050[[ii]] <- omat9
	out.ours075[[ii]] <- omat10
	out.ours095[[ii]] <- omat11
	out.ours100[[ii]] <- omat12
}

caron005.mn <- Reduce('+', out.caron005)/nsim
caron025.mn <- Reduce('+', out.caron025)/nsim
caron050.mn <- Reduce('+', out.caron050)/nsim
caron075.mn <- Reduce('+', out.caron075)/nsim
caron095.mn <- Reduce('+', out.caron095)/nsim
caron100.mn <- Reduce('+', out.caron100)/nsim

ours005.mn <- Reduce('+', out.ours005)/nsim
ours025.mn <- Reduce('+', out.ours025)/nsim
ours050.mn <- Reduce('+', out.ours050)/nsim
ours075.mn <- Reduce('+', out.ours075)/nsim
ours095.mn <- Reduce('+', out.ours095)/nsim
ours100.mn <- Reduce('+', out.ours100)/nsim

zlim <- range(c(c(caron005.mn),c(caron025.mn),c(caron050.mn),c(caron075.mn),c(caron095.mn),c(caron100.mn),
                c(ours005.mn),c(ours025.mn),c(ours050.mn),c(ours075.mn),c(ours095.mn),c(ours100.mn)))


# Figure 1 in the main document
# pdf("~/Research/BYU/SpaceTimePPM/latex/plots/caron3.pdf", height=7, width=10)
par(mfrow=c(2,3))
image.plot(caron005.mn, zlim=zlim, axes=FALSE, main=expression(alpha==0.05))
mtext(text=c(paste("",1:nt)), side=2, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)
mtext(text=c(paste("",1:nt)), side=1, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)

image.plot(caron025.mn, zlim=zlim, axes=FALSE, main=expression(alpha==0.25))
mtext(text=c(paste("",1:nt)), side=2, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)
mtext(text=c(paste("",1:nt)), side=1, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)

image.plot(caron050.mn, zlim=zlim, axes=FALSE, main=expression(alpha==0.50))
mtext(text=c(paste("",1:nt)), side=2, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)
mtext(text=c(paste("",1:nt)), side=1, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)

image.plot(caron075.mn, zlim=zlim, axes=FALSE, main=expression(alpha==0.75))
mtext(text=c(paste("",1:nt)), side=2, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)
mtext(text=c(paste("",1:nt)), side=1, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)

image.plot(caron095.mn, zlim=zlim, axes=FALSE, main=expression(alpha==0.95))
mtext(text=c(paste("",1:nt)), side=2, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)
mtext(text=c(paste("",1:nt)), side=1, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)

image.plot(caron100.mn, zlim=zlim, axes=FALSE, main=expression(alpha==1.0))
mtext(text=c(paste("",1:nt)), side=2, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)
mtext(text=c(paste("",1:nt)), side=1, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)
# dev.off()


# Figure 2 in the main document
# pdf("~/Research/BYU/SpaceTimePPM/latex/plots/ours3.pdf", height=7, width=10)
par(mfrow=c(2,3))
image.plot(ours005.mn, zlim=zlim,  axes=FALSE, main=expression(alpha==0.05))
mtext(text=c(paste("",1:nt)), side=2, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)
mtext(text=c(paste("",1:nt)), side=1, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)

image.plot(ours025.mn, zlim=zlim, axes=FALSE, main=expression(alpha==0.25))
mtext(text=c(paste("",1:nt)), side=2, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)
mtext(text=c(paste("",1:nt)), side=1, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)

image.plot(ours050.mn, zlim=zlim, axes=FALSE, main=expression(alpha==0.50))
mtext(text=c(paste("",1:nt)), side=2, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)
mtext(text=c(paste("",1:nt)), side=1, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)

image.plot(ours075.mn, zlim=zlim, axes=FALSE, main=expression(alpha==0.75))
mtext(text=c(paste("",1:nt)), side=2, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)
mtext(text=c(paste("",1:nt)), side=1, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)

image.plot(ours095.mn, zlim=zlim, axes=FALSE, main=expression(alpha==0.95))
mtext(text=c(paste("",1:nt)), side=2, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)
mtext(text=c(paste("",1:nt)), side=1, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)

image.plot(ours100.mn, zlim=zlim, axes=FALSE, main=expression(alpha==1))
mtext(text=c(paste("",1:nt)), side=2, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)
mtext(text=c(paste("",1:nt)), side=1, line=0.3, at=seq(0,1,length=nt), las=1, cex=0.8)
# dev.off()



# Function that cares out Monte Carlo experiment that produces Figure 3 in the main document
MCsim <- function(N, Tm, M, nRep, tau, sig, phi0=0, phi1=0, rhos,TYPE="Random"){


	corMatP <- matrix(NA, N, Tm)
	corMatPmn <- matrix(NA, length(rhos), Tm)
	corMatPsd <- matrix(NA, length(rhos), Tm)
	Ylist <- list()

	theta <- 0

	for(ii in 1:length(rhos)){
	

		cat("ii = ", ii, "\n")
		print(date())
		scorMatP <- matrix(0,Tm, Tm) 
	
		for(jj in 1:nRep){

			Pout <- rtpartition1(N=N,M=M,rho=rhos[ii],ntime=Tm,tau=tau,sig=sig,Caron=FALSE,FirstPart=NULL,TYPE=TYPE,
			                    phi0=phi0, phi1=phi1)
		
			Ylist[[jj]] <- Pout$YMat

		}

		for(t in 1:Tm){
			for(k in 2:N){
				tmp <- t(sapply(Ylist, function(x) x[c(1,t),k])) 
				corMatP[k,t] <- cor(tmp[,1], tmp[,2])
			}
		}
		corMatPmn[ii,] <- apply(corMatP[-1,],2,mean)
		corMatPsd[ii,] <- apply(corMatP[-1,],2,sd)

	}

    # Note that here I am writing out to file of the working directory
	write.table(corMatPmn, paste("corMatmn",TYPE,phi1,".txt", sep=""), col.name=1:Tm, row.name=rhos)
	write.table(corMatPsd, paste("corMatsd",TYPE,phi1,".txt", sep=""), col.name=1:Tm, row.name=rhos)
}

set.seed(1)



rhos <- c(0,0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
nRep <- 10000
MCsim(N=25, Tm=10, M=1, nRep=nRep, tau=1, sig=0.1, phi1=0, rhos=rhos, TYPE="Random")
MCsim(N=25, Tm=10, M=1, nRep=nRep, tau=1, sig=0.1, phi1=0, rhos=rhos, TYPE="Deterministic")
MCsim(N=25, Tm=10, M=1, nRep=nRep, tau=1, sig=0.1, phi1=0.01, rhos=rhos, TYPE="AR1")
MCsim(N=25, Tm=10, M=1, nRep=nRep, tau=1, sig=0.1, phi1=0.1, rhos=rhos, TYPE="AR1")
MCsim(N=25, Tm=10, M=1, nRep=nRep, tau=1, sig=0.1, phi1=0.25, rhos=rhos, TYPE="AR1")
MCsim(N=25, Tm=10, M=1, nRep=nRep, tau=1, sig=0.1, phi1=0.5, rhos=rhos, TYPE="AR1")
MCsim(N=25, Tm=10, M=1, nRep=nRep, tau=1, sig=0.1, phi1=0.75, rhos=rhos, TYPE="AR1")
MCsim(N=25, Tm=10, M=1, nRep=nRep, tau=1, sig=0.1, phi1=0.9, rhos=rhos, TYPE="AR1")
MCsim(N=25, Tm=10, M=1, nRep=nRep, tau=1, sig=0.1, phi1=0.99, rhos=rhos, TYPE="AR1")


# Read in results produce from the above code
# Repalce "marginal_correlations" to folder where output from above function calls is stored
files <- list.files("marginal_correlations")
corMatsmn <- list()
corMatssd <- list()
for(ii in 1:(length(files)/2)){
  corMatsmn[[ii]] <- read.table(paste("/marginal_correlations/", files[ii],sep=""),
                           col.names=paste("lag", 0:9, sep=""))
  corMatssd[[ii]] <- read.table(paste("/marginal_correlations/",files[ii+(length(files)/2)],sep=""),
                           col.names=paste("lag", 0:9, sep=""))
}


phi <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1, 0)
rhos <- c(0,0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)


# Code that produces Figure 3 of the main document
# pdf("~/Research/BYU/SpaceTimePPM/latex/plots/MarginalCorrelationPlots5.pdf", height=8, width=13)

  par(mfrow=c(2,3))

  for(ii in c(9, 2, 4, 6, 7, 8)){

    plot(t(corMatsmn[[ii]][1,]), ylim=c(0.0,1), ylab="auto correlation", xlab="lag", 
			main=bquote(phi[1] == .(phi[ii])), type='n', cex.lab=1.25, cex.axis=1.25, cex.main=1.5)

    jj <- 1
    for(k in c(1, 4, 5, 6,7)){
      lines(t(corMatsmn[[ii]][k,]), type="b", col="black", lty=jj, lwd=2)
      jj <- jj + 1
    }

    legend(x="topright", 
        legend=c(as.expression(bquote(alpha == .(rhos[1]))), 
		                       bquote(alpha == .(rhos[4])),
	    	                   bquote(alpha == .(rhos[5])),
	        	               bquote(alpha == .(rhos[6])),
	            	           bquote(alpha == .(rhos[7])),
	            	           bquote(alpha == .(rhos[8]))), 
        ncol=2, col="black", lty=1:5, pch=1, lwd=2,
        cex=1.5, seg.len=2.5)
  }
# dev.off()

