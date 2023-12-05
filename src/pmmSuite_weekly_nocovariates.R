library(ppmSuite)
load("../data/df_weekly.Rdata")

nobs <- 53 
nsubject <- length(unique(df_weekly$IDStations))

df_vuoto <- data.frame(matrix(nrow = nsubject, ncol = nobs))


y = matrix(nrow = 0, ncol = nobs)
for (i in unique(df_weekly$IDStations))
{
   y <- rbind(y, filter(df_weekly, IDStations == i)$AQ_pm10)
  
}
 
dat <- data.frame(y=c(y),
                   z=rep(1:nobs, times=nsubject),
                   Name=rep(1:nsubject, each=nobs))

subject_obs_vec <- dat$Name
nknots <- 15 
# Small number of iterates for illustrative purposes only
niter <- 5000
nburn <- 2000
nthin <- 3
nout <- (niter-nburn)/nthin
z <- dat$z
## the order here is c(mu0, s20, v, k0, nu0, a0, alpha)
## If simularity is N-NIG then k0 and nu0 are used but v is not
## If simularity is N-N then v is used but no k0 and nu0
simparms <- c(0.0, 1.0, 0.1, 1.0, 1.0, 0.1, 1)
fits <- list()
# fit vgrf only
y <- dat$y
modelPriors <- c(0.5, # Asig
                 1000^2, # s2_mu
                 0, # mb0
                 1000^2, # s2b0
                 1, # as2b0
                 1, # bs2b0
                 1, # at
                 1.0/0.05) # bt

fit <- curve_ppmx(y=cbind(y), z=z,
                  subject=subject_obs_vec,
                  Xcon = NULL, Xcat = NULL,
                  Xconp=NULL, Xcatp=NULL,
                  PPM=TRUE, M=1,
                  q=3, rw_order=1, balanced=1,
                  nknots=nknots, npredobs=1,
                  Aparm=100,
                  modelPriors=modelPriors,
                  similarity_function=1,
                  consim=1, calibrate=0,
                  simParms=simparms,
                  mh=c(0.1, 1e-4),
                  draws=niter,
                  burn=nburn,
                  thin=nthin)
Hmat <- fit$Hmat
# For a point estimate of partition, take first MCMC interate
# This is done only for illustrative purposes. Recommend using
# the salso R package.
p.est <- fit$Si[1,]
nc <- length(unique(p.est))
oldpar <- par(no.readonly = TRUE)
# Plot individual subject fits.
tmp <- c(1,6,11,16)
windows()
par(mfrow=c(2,2))
for(j in tmp){
  bmn <- apply(fit$beta[j,,],1,mean)
  b0mn <- mean(fit$beta0[,j])
  ytmp <- y[dat$Name==j]
  b0vec <- rep(b0mn, nobs)
  plot(1:nobs,c(ytmp),
       type='n',ylab="Response",
       xlab="Time")
  points(1:nobs,ytmp)
  lines(1:nobs, b0vec+Hmat%*%bmn, col=p.est[j],lwd=2)
}
# plot all curves in one plot
windows()
par(mfrow=c(1,1))
plot(dat$z, dat$y, type="n",ylab="",xlab="Time")
for(j in 1:nsubject){
  bmn <- apply(fit$beta[j,,],1,mean)
  b0mn <- mean(fit$beta0[,j])
  b0vec <- rep(b0mn, nobs)
  lines((1:nobs), b0vec+Hmat%*%bmn, col=p.est[j],lwd=0.5)
}
par(oldpar)

salso(t(fit$Si),structure="clustering",loss="binder")

