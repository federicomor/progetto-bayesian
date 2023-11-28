# Thurs. April. 15 2021
# Rerun this simulation with updated salso and drpm packages
#
# Try and generate data from the model and then recover it
# by fitting synthetic data to model.

# To generate data from the independent cluster-specific parameters model
# the rtpartition function created in Functions.R file


source("Functions.R")
library(drpm)
library(salso)
library(mclust)
library(TeachingDemos)

N <- 50; Tm<-5; M<-1;
alpha <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.9999)
#alpha <- c(0.5)

# simulation study with 2 data sets
ndata <- 2
nalpha <- length(alpha)


for(jj in 1:nalpha){

  cat("alpha = ", alpha[jj], "\n")
  set.seed(1 + 100*alpha[jj])
  cat("set.seed = ", 1 + 100*alpha[jj], "\n")

  cat("jj = ", jj, "\n")

  out.cov <- NULL
  out.rand <- matrix(NA, nrow=ndata, ncol=Tm)
  out.mean.cov <- matrix(NA, nrow=ndata, ncol=1)
  out.alpha.cov <- matrix(NA,nrow=ndata, ncol=1)
  out.adjust.rand.first.last <- matrix(NA, nrow=ndata, ncol=1)
  out.adjust.rand.cont <- matrix(NA, nrow=ndata, ncol=Tm-1)

  fit.out <- matrix(NA, nrow=ndata, ncol=4)
  colnames(fit.out)<-c("DepPartLPML","DepPartWAIC","IndPartLPML","IndPartWAIC")

  for(ii in 1:ndata){
    cat("ii = ", ii, "\n")

    dat <- rtpartition1(N=N,M=M,rho=alpha[jj],ntime=Tm,tau=5,sig=1,
						Caron=FALSE,FirstPart=NULL,TYPE="Random",
					    phi0=0, phi1=1)

    y <- t(dat$YMat)
    ciMat <- dat$ciMat
    muMat <- dat$muMat


    niter=10000; nburn=5000; nthin=5


    #             # m0, s20, A, At,  Al,  at, bt, be
    modelPriors <- c(0, 100, 5, 10,  10,  1,  1,  1)

    # Dependent CRP
    out <- drpm_fit(y=y, 
    				#global_alpha=TRUE,
                    alpha_0 = FALSE,
                    eta1_0 = TRUE,
                    phi1_0 = TRUE,
                    modelPriors=modelPriors,
                    draws=niter, burn=nburn, thin=nthin,verbose=TRUE)

    fit.out[ii,1] <- out$lpml
    fit.out[ii,2] <- out$waic

    mu.int <-  aperm(apply(out$mu,c(1,2),emp.hpd), c(3,1,2))

    rho <- list()

    for(k in 1:Tm){
      rho[[k]] <- salso(t(out$Si[k,,]), loss="binder")
      out.rand[ii,k] <- adjustedRandIndex(rho[[k]], ciMat[k,])

      out.cov[k] <- sum(apply(mu.int[,,k] - muMat[k,],1,prod)	< 0)
    }

    out.adjust.rand.first.last[ii,] <- adjustedRandIndex(rho[[1]], rho[[k]])

    for(k in 1:(Tm-1)){
      out.adjust.rand.cont[ii,k] <- adjustedRandIndex(rho[[k]], rho[[k+1]])
    }

    out.mean.cov[ii] <- sum(out.cov)/(N*Tm)

    out.alpha.cov[ii] <- prod(emp.hpd(out$alpha,0.95) - alpha[jj]) < 0

    # Indepdent CRP
    out1 <- drpm_fit(y=y, 
    				#global_alpha=TRUE,
                    alpha_0 = TRUE,
                    eta1_0 = TRUE,
                    phi1_0 = TRUE,
                    modelPriors=modelPriors,
                    draws=niter, burn=nburn, thin=nthin,verbose=TRUE)

    fit.out[ii,3] <- out1$lpml
    fit.out[ii,4] <- out1$waic

    write.table(fit.out, paste("fit_", alpha[jj], ".txt",sep=""), row.names=FALSE, col.names=TRUE)

    # I write to file the individual model fits, but this is disk-space consuming.
#    save.image(paste("simstudy1fits/alpha_", alpha[jj], "_data_", ii,".RData",sep=""))

    write.table(round(out.rand,4),
      paste("simstudy1results/out.rand_", alpha[jj],".txt", sep=""), row.names=FALSE, col.names=FALSE)

    write.table(round(out.mean.cov,4),
      paste("simstudy1results/out.mean.cov_", alpha[jj], ".txt",sep=""), row.names=FALSE, col.names=FALSE)

    write.table(round(out.alpha.cov,4),
      paste("simstudy1results/out.alpha.cov_", alpha[jj], ".txt",sep=""), row.names=FALSE, col.names=FALSE)

    write.table(round(out.adjust.rand.first.last,4),
      paste("simstudy1results/out.rand.first.last_", alpha[jj], ".txt",sep=""), row.names=FALSE, col.names=FALSE)

    write.table(round(out.adjust.rand.cont,4),
      paste("simstudy1results/out.rand.cont_", alpha[jj], ".txt",sep=""), row.names=FALSE, col.names=FALSE)
	
  }

}



# Organize simulation 1 results by creating table 2

# This should be set to directory where output is stored. 
# "simstudy1results" is provided with output from our run. 
# setwd("simstudy1results")

files <- list.files("simstudy1results")

out.mn <- list()
out.sd <- list()
out.dat <- list()
fit.lpml <- list()
fit.waic <- list()

for(i in 1:length(files)){
  print(files[i])
  if(i < 8){
    tmp <- read.table(files[i], header=TRUE)
    fit.lpml[[i]] <- apply(tmp,2,function(x) mean(x[is.finite(x)]))[c(1,3)]
    fit.waic[[i]] <- apply(tmp,2,function(x) mean(x[is.finite(x)]))[c(2,4)]
  } else {
    tmp <- read.table(files[i], header=FALSE)
    out.mn[[i-7]] <- apply(cbind(tmp),2,mean, na.rm=TRUE)
    out.sd[[i-7]] <- apply(cbind(tmp),2,function(x) {n <- sum(!is.na(x)); sd(x[!is.na(x)]/sqrt(n))})
    out.dat[[i-7]] <- tmp
  }
}
ord <- c(7,1:6)
covs.mn <- cbind(out.mn[ord],out.sd[ord], out.mn[ord+7],out.sd[ord+7]); 
colnames(covs.mn) <- c("alpha","alphase","y","yse")
rand.true <- cbind(matrix(unlist(out.mn[ord+14]), nrow=length(alpha), byrow=TRUE),
                   matrix(unlist(out.sd[ord+14]), nrow=length(alpha), byrow=TRUE))
rand.cont <- cbind(matrix(unlist(out.mn[ord+21]), nrow=length(alpha), byrow=TRUE), out.mn[ord+28], 
                   matrix(unlist(out.sd[ord+21]), nrow=length(alpha), byrow=TRUE), out.sd[ord+28])
lpml <- matrix(unlist(fit.lpml[ord]), nrow=length(alpha), byrow=TRUE)
waic <- matrix(unlist(fit.waic[ord]), nrow=length(alpha), byrow=TRUE)

library(xtable)

# Table 2 of the paper
xtable(cbind(rand.cont[,c(1,6,5,10)],covs.mn, round(waic)), digits=2)


# Table 3 in the supplementary material
xtable(rand.true[,c(1,6,2,7,3,8,4,9,5,10)])

# plot adjusted rand for \rho_1 and \rho_2 as a function of alpha
# Not included in final draft of paper.
val <- c(out.dat[[28]][,1], out.dat[[22]][,1], out.dat[[23]][,1],
         out.dat[[24]][,1], out.dat[[25]][,1], out.dat[[26]][,1],
         out.dat[[27]][,1])
grp <- rep(c(0,0.1,0.25,0.5,0.75,0.9,0.9999), each=100)
x11()
boxplot(val~grp, ylab="", xlab=expression(alpha))
mtext(expression("Adjusted Rand Index between"~hat(rho)[1]~"and"~hat(rho)[2]), side=2, line=2.5)
