#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

# To excute code, from command line type the following
# nohup ./simstudy3_ref.R 100 4 100 5 0.0 0.5 > out85.txt &

# This simulation study relies heavily on code that I've
# included in the drpm R package.  These fit the linearDDP
# weighted DDP, and DP models.

# args <- c(100, 4, 100, 5, 0.0, 0.5)

ndata <- as.numeric(args[1])
datatype <- as.numeric(args[2])
N <- as.numeric(args[3])
Tm <- as.numeric(args[4])
correlation <- as.numeric(args[5])
stdev <- as.numeric(args[6])

cat("ndata = ", ndata, "\n")
cat("N = ", N, "\n")
cat("datatype = ", datatype, "\n")
cat("Tm = ", Tm, "\n")
cat("correlation = ", correlation, "\n")
cat("stdev = ", stdev, "\n")

niter <- 50000
nburn <- 10000
nthin <- 20
nout <- (niter - nburn)/nthin

library(drpm)
library(BNPmix) #DDPdensity
library(splines) #spline.des
library(salso)
library(mclust)

# This simulation study generates data from a AR(1) process
# using 5 or 10 time points and 100 observations at each one.
# Then to each data set, we fit the following methods.
#
# 1. our method based on model 5 of our paper
#    (i.e., temporal dependence in partition model only)
#
# 2. weighted dependent dirichlet process code created by me
#
# 3. linear dependent dirichlet process code created by me
#
# 4. CRP independent at each time point code created by me
#
# 5. CRP based on nT observations code created by me.
#
# 6. The univariate Griffiths-Milne dependent Dirichlet
#    process mixture model with Gaussian kernel,
#    for partially exchangeable data This is
#    fit using the BNPmix package

# Function to generate data
time.dat.gen <- function(N, Tm, auto_cor=0.5, sd=1){
  if(auto_cor == 0){
        YMat <- t(replicate(N, arima.sim(model=list(), n=Tm, sd=sd)))
  } else {
    YMat <- t(replicate(N, arima.sim(model=list(ar=c(auto_cor)), n=Tm, sd=sd)))
  }
  YMat
}

# Function that creates B-spline basis
bspline <-function(X., XL., XR., NDX., BDEG.){

	dx <- (XR. - XL.)/NDX.
	knots <- seq(XL. - BDEG.*dx, XR. + BDEG.*dx, by=dx)
	B <- spline.des(knots, X., BDEG.+1, 0*X.)$design
	res <- list(B = B, knots = knots)
	res
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


# Function that relabels component indicators with using
# the Blocked GIBBS sampler of weighted DDP and DP models
relabel = function(z){
  labs = rep(1:length(unique(z)),times=table(z))
  z[order(z)] = labs
  z
}


waic <- matrix(NA, nrow=ndata, ncol=6)
ari <- matrix(NA, nrow=ndata, ncol=6)
lpml <- matrix(NA, nrow=ndata, ncol=6)
colnames(waic) <- colnames(ari) <- colnames(lpml) <- c("drpm", "wddp", "lddp", "gmddp", "indcrp","longcrp")
probpeg <- matrix(NA, nrow=ndata,ncol=2)
colnames(probpeg) <- c("gamma","alpha")

# ii <- 1
for(ii in 1:ndata){

  cat("dataset ================================================= ", ii, "\n")
  set.seed(100 + ii);
  cat("seed = ", 100+ii, "\n")

  # This data type studies performance when no underlying clustering exists
  # but there is autocorrelation in the observations
  if(datatype == 1){
    mn <- c(0, 0, 0, 0)
    cor <- rep(correlation, 4)
    Ymat <- matrix(NA, nrow=N, ncol=Tm)
    for(c in 1:4){
      Ymat[((c-1)*(N/4) +1):(c*(N/4)),1:Tm] <- mn[c] + time.dat.gen(N=N/4, Tm=Tm, auto_cor=cor[c], sd=stdev)
    }
    tp <- matrix(rep(1:4, each=(N/4)*Tm), ncol=Tm, byrow=TRUE)
  }

  # This datatype has clusters based on mean and clusters composition (i.e., who are clustered
  # together) changes over time.  Also have autocorrelation.
  if(datatype == 2){
    mn <- c(-2, 0, 2, 4)
    cor <- rep(correlation, 4)
    Ymat <- NULL
    tp <- cbind(rep(1:4, each=N/4))
    for(t in 1:Tm){
      Ymat <- cbind(Ymat, t(mn[tp[,t]] + time.dat.gen(N=N, Tm=1, auto_cor=correlation, sd=stdev)))
      if(t == Tm) break
      tp <- cbind(tp, c(tp[,t][-c(1:4)], tp[,t][c(1:4)]))
    }
  }

  # This datatype is like datatype 2 but clusters mimic DP rich get richer type.
  if(datatype == 3){
    mn <- c(-2, 0, 2, 4)
    cor <- rep(correlation, 4)
    csize <- c(55, 25, 15, 5)
    Ymat <- NULL
    tp <- cbind(rep(1:4, times=csize))
    for(t in 1:Tm){
      Ymat <- cbind(Ymat, t(mn[tp[,t]] + time.dat.gen(N=N, Tm=1, auto_cor=correlation, sd=stdev)))
      if(t == Tm) break
      tp <- cbind(tp, c(tp[,t][-c(1:4)], tp[,t][c(1:4)]))
    }
  }

  # This datatype has clusters but based on autocorrelation (our model isn't equiped to detect this).
  if(datatype == 4){
    mn <- c(0, 0, 0, 0)
    cor <- c(-0.75, -0.25, 0.25, 0.75)
    Ymat <- matrix(NA, nrow=N, ncol=Tm)
    for(c in 1:4){
      Ymat[((c-1)*(N/4) +1):(c*(N/4)),1:Tm] <- mn[c] + time.dat.gen(N=N/4, Tm=Tm, auto_cor=cor[c], sd=stdev)
    }
    tp <- matrix(rep(1:4, each=(N/4)*Tm), ncol=Tm, byrow=TRUE)
  }


  Tmat <- matrix(rep(1:Tm, each=N), nrow=N, byrow=FALSE)

  Yvec <- 1*c(t(Ymat))
  Tvec <- cbind(c(t(Tmat)))
  part.long <- c(t(tp))


  ##################################################################
  #
  # our method
  #             # m0, s20, A,            At,   Al,   at, bt, be
  modelPriors <- c(0, 100, 0.5*sd(Yvec), 100,  100,  1,  1,  1)

  drpm1 <- drpm_fit(y=Ymat, global_alpha=FALSE,
                 alpha_0 = FALSE,
                 eta1_0 = TRUE,
                 phi1_0 = TRUE,
                 modelPriors=modelPriors,
                 draws=niter, burn=nburn, thin=nthin,verbose=TRUE)

  probpeg[ii,1] <-  mean(apply(apply(drpm1$gamma, c(1,2), mean),1,mean)[-1])
  probpeg[ii,2] <-  mean(apply(drpm1$alpha,2,mean)[-1])

  ari.mat.drpm <- matrix(NA, nout, Tm)
  for(t in 1:Tm){
    for(jj in 1:nout){
      ari.mat.drpm[jj,t] <- adjustedRandIndex(drpm1$Si[t,,jj], tp[,t])
    }
  }

  ##################################################################
  #
  # Weighted DDP  Note we are adding jitter to T to improve performace
  #           #  m  s2   nu1 Psi - rate
  priorparms <- c(0, 25,  4,  1.0);
  wddp <- weightedDDP(draws=niter, burn = nburn, thin=nthin,
                      y=Yvec, x=jitter(Tvec, amount=0.5), ygrid=1,
                      xpred=1, priorparms=priorparms, N=30, Mdp=1)

  ari.mat.wddp <- matrix(NA, nout, Tm)
  for(t in 1:Tm){
    for(jj in 1:nout){
      ari.mat.wddp[jj,t] <- adjustedRandIndex(wddp$Si[jj, (N*(t-1)+1):(N*t)], tp[,t])
    }
  }



  ##################################################################
  #
  # Linear DDP
  #
  Tbs <- cbind(bspline(seq(min(Tvec)-0.05,max(Tvec)+0.05,
                           length=length(Yvec)),
                       min(Tvec[,1])-0.1,max(Tvec[,1])+0.1,10,3)$B)
        # mb0  s2b0  scaleSig0, df Sig0
  pp <- c(0,  25,    10,        dim(Tbs)[2]+2, 1)

  lddp <- linearDDP(draws=niter, burn=nburn, thin=nthin, y=Yvec, x=Tbs,
                      ygrid=NULL, xpred=NULL, priorparms=pp, Mdp=1)

  ari.mat.lddp <- matrix(NA, nout, Tm)
  for(t in 1:Tm){
    for(jj in 1:nout){
      ari.mat.lddp[jj,t] <- adjustedRandIndex(lddp$Si[jj, (N*(t-1)+1):(N*t)], tp[,t])
    }
  }



  ##################################################################
  #
  # GM - DDP
  #
  grid <- Yvec
  gmddp <- DDPdensity(y = Yvec, group = c(Tvec),
                prior = list(strength=1, wei=0.5, m0=mean(Yvec),k0=1, a0=2, b0=0.5*sd(Yvec)),
                mcmc = list(niter = niter, nburn = nburn),
                output=list(grid=grid, out_type="FULL"))
  llike <- NULL
  thin <- seq(1, niter-nburn, by = nthin)
  for(t in 1:Tm){
    llike <- cbind(llike, t(log(gmddp$density[((t-1)*N + 1):(t*N), t,thin])))
  }
  gmddp_waic <- lpml.robust(llike)[6]
  gmddp_lpml <- lpml.robust(llike)[1]

  gmddp_Si <- gmddp$clust[thin,]
  ari.mat.gmddp <- matrix(NA, nout, Tm)
  for(t in 1:Tm){
    for(jj in 1:nout){
      ari.mat.gmddp[jj,t] <- adj.rand.index(gmddp_Si[jj, (N*(t-1)+1):(N*t)], tp[,t])
    }
  }



  ##################################################################
  #
  # For the two extremes associated with Caron's method,
  # I created another set of code that fits a DPM model
  # This means that the theta parameter for in our method
  # is set to a fixed value (0). I also use an IG prior on
  # the component variance.
  #
  # independent CRPs
  #
  #
  indcrpllike <- NULL
  waic.tmp <- NULL
  ari.mat.indcrp <- matrix(NA, nrow=nout, ncol=Tm)
  indcrpfits <- list()
  for(t in 1:Tm){

    indcrp <- gaussian_crp(Ymat[,t], m=0, v=10^2, A=0.5*sd(Yvec), A0=100,
                          mh=c(0.5, 0.5), alpha=1,
                          niter=niter, nburn=nburn, nthin=nthin)
    indcrpfits[[t]] <- indcrp
    indcrpllike <- cbind(indcrpllike, indcrp$llike)
    waic.tmp <- c(waic.tmp, indcrp$waic)
    for(jj in 1:nout){
      ari.mat.indcrp[jj,t] <- adjustedRandIndex(indcrp$Si[jj,], tp[,t])
    }
  }
  indcrp_waic <- lpml.robust(indcrpllike)[6]
  indcrp_lpml <- lpml.robust(indcrpllike)[1]


  ##################################################################
  #
  # one long CRP
  #
  # This code fits the one CRP model using our function
  # but it is quite slow.  So I run a DPM model using the
  # concatinated response vector.
  #


  longcrp <- gaussian_crp(Yvec, m=0, v=10^2, A=0.5*sd(Yvec), A0=100,
                         mh=c(0.5, 0.5), alpha=1,
                         niter=niter, nburn=nburn, nthin=nthin)

  pe.longcrp <- salso(t(apply(longcrp$Si,1,relabel)), loss="VI", nCores=1)

  ari.mat.longcrp <- matrix(NA, nout, Tm)
  for(t in 1:Tm){
    for(jj in 1:nout){
      ari.mat.longcrp[jj,t] <- adjustedRandIndex(longcrp$Si[jj, (N*(t-1)+1):(N*t)], tp[,t])
    }
  }


  waic[ii,1] <- drpm1$waic
  waic[ii,2] <- wddp$waic
  waic[ii,3] <- lddp$waic
  waic[ii,4] <- gmddp_waic
  waic[ii,5] <- indcrp_waic
  waic[ii,6] <- longcrp$waic

  lpml[ii,1] <- drpm1$lpml
  lpml[ii,2] <- wddp$lpml
  lpml[ii,3] <- lddp$lpml
  lpml[ii,4] <- gmddp_lpml
  lpml[ii,5] <- indcrp_lpml
  lpml[ii,6] <- longcrp$lpml

  ari[ii,1] <- mean(apply(ari.mat.drpm,2,mean))
  ari[ii,2] <- mean(apply(ari.mat.wddp,2,mean))
  ari[ii,3] <- mean(apply(ari.mat.lddp,2,mean))
  ari[ii,4] <- mean(apply(ari.mat.gmddp,2,mean))
  ari[ii,5] <- mean(apply(ari.mat.indcrp,2,mean))
  ari[ii,6] <- mean(apply(ari.mat.longcrp,2,mean))

   # This is for running it locally

  # This is for running on the server
  dir <- paste0("simstudy3results/datatype_",
                      datatype, "_corr_",correlation, "_stdev_", stdev,
                     "_time_", Tm)

  write.table(waic, paste0(dir, "_waic.txt"),
                col.names=TRUE, row.names=FALSE)

  write.table(lpml, paste0(dir, "_lpml.txt"),
                col.names=TRUE, row.names=FALSE)

  write.table(ari, paste0(dir, "_ari.txt"),
                col.names=TRUE, row.names=FALSE)

  write.table(probpeg, paste0(dir,  "_probpeg.txt"),
                col.names=TRUE, row.names=FALSE)


  rm(list=c("longcrp", "indcrpllike", "indcrp", "lddp",
            "llike", "gmddp","wddp","drpm1", "gmddp_waic",
            "grid","indcrp_waic","jj","modelPriors","pp","t","Tbs",
            "pe.drpm","pe.gmddp","pe.indcrp","pe.lddp","pe.longcrp",
            "pe.wddp","priorparms","thin","waic.tmp", "indcrp_lpml",
            "ari.mat.longcrp","ari.mat.indcrp","ari.mat.drpm",
            "ari.mat.wddp","ari.mat.lddp","ari.mat.gmddp"))



}



# The code that follows reads in results from simulation study, organizes it,
# and then produces tables and plots that are included in the main article of
# the paper.
#
# Note that the "dir" is the directory where simulation output is stored.

  # This code organizes simulation output that was stored in the folder "simstudy3results"
  dir <- "~/directory to folder where output is saved/simstudy3results/"

  files <- list.files(dir)


  ari.mat <- data.frame(corr=NULL, stdev=NULL, datatype=NULL, time=NULL, procedure=NULL, ari=NULL)
  waic.mat <- data.frame(corr=NULL, stdev=NULL, datatype=NULL, time=NULL, procedure=NULL, waic=NULL)
  lpml.mat <- data.frame(corr=NULL, stdev=NULL, datatype=NULL, time=NULL, procedure=NULL, waic=NULL)
  prob.peg <- data.frame(corr=NULL, stdev=NULL, datatype=NULL, time=NULL, probpeg=NULL)
  for(jj in 1:length(files)){
    factors <- strsplit(files[jj], "\\_")[[1]][c(2,4,6,8)]
    datatype <- as.numeric(factors[1])
    corr <- as.numeric(factors[2])
    stdev <- as.numeric(factors[3])
    time <- as.numeric(factors[4])
    metric <- strsplit(files[jj], "\\_|\\.")[[1]]
    metric <- metric[length(metric) - 1]

    tmp <- read.table(paste0(dir, files[jj]), header=TRUE)

    nproc <- ncol(tmp)
    ndata <- nrow(tmp)
    if(metric == "ari"){
      ari.mat <- rbind(ari.mat,
                       data.frame(corr=rep(corr, nproc*ndata),
                                  stdev = rep(stdev, nproc*ndata),
                                  datatype = rep(datatype, nproc*ndata),
                                  time = rep(time, nproc*ndata),
                                  procedure = rep(colnames(tmp), each=ndata),
                                  ari = c(unlist(tmp))))
    }
    if(metric == "waic"){
      waic.mat <- rbind(waic.mat,
                        data.frame(corr=rep(corr, nproc*ndata),
                                   stdev = rep(stdev, nproc*ndata),
                                   datatype = rep(datatype, nproc*ndata),
                                   time = rep(time, nproc*ndata),
                                   procedure = rep(colnames(tmp), each=ndata),
                                   waic = c(unlist(tmp))))

    }
    if(metric == "lpml"){
      lpml.mat <- rbind(lpml.mat,
                        data.frame(corr=rep(corr, nproc*ndata),
                                   stdev = rep(stdev, nproc*ndata),
                                   datatype = rep(datatype, nproc*ndata),
                                   time = rep(time, nproc*ndata),
                                   procedure = rep(colnames(tmp), each=ndata),
                                   waic = c(unlist(tmp))))

    }
    if(metric == "probpeg"){
      prob.peg <- rbind(prob.peg,
                        data.frame(corr=rep(corr, ndata),
                                  stdev = rep(stdev,ndata),
                                  datatype = rep(datatype, nproc*ndata),
                                  time = rep(time, nproc*ndata),
                                  probpeg = tmp[,1]))

    }

  }

  row.names(waic.mat) <- row.names(ari.mat) <- row.names(lpml.mat) <- row.names(prob.peg) <- NULL
  waic.mat$procedure <- factor(waic.mat$procedure,
                               levels = c("drpm", "indcrp", "longcrp", "lddp", "wddp", "gmddp"),
                               labels = c("drpm", "ind_crp", "static_crp", "lddp", "wddp", "gmddp"))
  waic.mat$time <- factor(waic.mat$time,
                               levels = c("5", "10"),
                               labels = c("T = 5", "T = 10"))
  waic.mat$stdev <- factor(waic.mat$stdev,
                               levels = c("0.5", "1"),
                               labels = c("v = 0.5", "v = 1"))
  lpml.mat$procedure <- factor(lpml.mat$procedure,
                               levels = c("drpm", "indcrp", "longcrp", "lddp", "wddp", "gmddp"),
                               labels = c("drpm", "ind_crp", "static_crp", "lddp", "wddp", "gmddp"))
  lpml.mat$time <- factor(lpml.mat$time,
                               levels = c("5", "10"),
                               labels = c("T = 5", "T = 10"))
  lpml.mat$stdev <- factor(lpml.mat$stdev,
                               levels = c("0.5", "1"),
                               labels = c("v = 0.5", "v = 1"))
  ari.mat$procedure <- factor(ari.mat$procedure,
                               levels = c("drpm", "indcrp", "longcrp", "lddp", "wddp", "gmddp"),
                               labels = c("drpm", "ind_crp", "static_crp", "lddp", "wddp", "gmddp"))
  ari.mat$time <- factor(ari.mat$time,
                               levels = c("5", "10"),
                               labels = c("T = 5", "T = 10"))
  ari.mat$stdev <- factor(ari.mat$stdev,
                               levels = c("0.5", "1"),
                               labels = c("v = 0.5", "v = 1"))
  data_mn <- function(data){
    metric <- strsplit(deparse(substitute(data)), "\\.")[[1]][1]
    data.mn <- aggregate(data[,6],
                  list(as.factor(data$corr), data$stdev,
                       as.factor(data$datatype), data$time,
                       data$procedure),
                median, na.rm=TRUE)
    colnames(data.mn) <- c("corr","stdev","datatype","time","procedure",metric)
    data.mn$corr <- as.numeric(as.character(data.mn$corr))
    data.mn$datatype <- as.numeric(as.character(data.mn$datatype))
    data.mn <- data.mn[data.mn$corr < 0.91,]
    data.mn
  }

  waic.mn <- data_mn(waic.mat)
  lpml.mn <- data_mn(lpml.mat)
  ari.mn <- data_mn(ari.mat)


  library(ggplot2)

  # datatype 1
  # Figure S.2 in the supplementary material
  gg.waic <- ggplot(data = waic.mn[waic.mn$datatype==1,], aes(x=as.numeric(corr), y=waic,col=procedure)) +
               geom_line() +
               geom_point() +
               facet_grid(time ~ stdev,scales="free") + theme_bw() +
               xlab("auto-correlation used to generate data") + ylab("WAIC")
#  ggsave("simstudy_datatype1_waic2.pdf", height=6, width=6)


  # datatype 2
  # Figure 4 in the main document
  gg.waic <- ggplot(data = waic.mn[waic.mn$datatype==2,], aes(x=as.numeric(corr), y=waic,col=procedure)) +
               geom_line() +
               geom_point() +
               facet_grid(time ~ stdev,scales="free") + theme_bw() +
          	  xlab("auto-correlation used to generate data") + ylab("WAIC")
#  ggsave("simstudy_datatype2_waic2.pdf", height=6, width=6)

  # Figure 4 in the main document
  gg.ari <- ggplot(data = ari.mn[ari.mn$datatype==2,], aes(x=as.numeric(corr), y=ari,col=procedure)) +
               geom_line() +
               geom_point() +
               facet_grid(time ~ stdev,scales="free") + theme_bw() +
          	  xlab("auto-correlation used to generate data") + ylab("Adjusted Rand Index")
#  ggsave("simstudy_datatype2_ari2.pdf", height=6, width=6)



  # datatype 3
  # These next two plots were not included in the document 
  # as the results are essentially identical to datatype 2
  gg.waic <- ggplot(data = waic.mn[waic.mn$datatype==3,], aes(x=as.numeric(corr), y=waic,col=procedure)) +
               geom_line() +
               geom_point() +
               facet_grid(time ~ stdev,scales="free") + theme_bw() +
          	  xlab("auto-correlation used to generate data") + ylab("WAIC")
#  ggsave("simstudy_datatype3_waic2.pdf", height=6, width=6)

  gg.ari <- ggplot(data = ari.mn[ari.mn$datatype==3,], aes(x=as.numeric(corr), y=ari,col=procedure)) +
               geom_line() +
               geom_point() +
               facet_grid(time ~ stdev,scales="free") + theme_bw() +
          	  xlab("auto-correlation used to generate data") + ylab("Adjusted Rand Index")
#  ggsave("simstudy_datatype3_ari2.pdf", height=6, width=6)


  # datatype 4
  # Figure S.3 in the supplmentary material
  gg.waic <- ggplot(data = waic.mat[waic.mat$datatype==4,], aes(x=procedure, y=waic)) +
               geom_boxplot() +
               facet_grid(time ~ stdev,scales="free") + theme_bw() + ylim(0, 4000) +
          	  xlab("procedure") + ylab("WAIC")
#  ggsave("simstudy_datatype4_waic2.pdf", height=6, width=8)

  # Figure S.4 in the supplmentary material
  gg.ari <- ggplot(data = ari.mat[ari.mat$datatype==4,], aes(x=procedure, y=ari)) +
               geom_boxplot() +
               facet_grid(time ~ stdev,scales="free") + theme_bw() +
          	  xlab("procedure") + ylab("Adjusted Rand Index")
#  ggsave("simstudy_datatype4_ari2.pdf", height=6, width=8)





