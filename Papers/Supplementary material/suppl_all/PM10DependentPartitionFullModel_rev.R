#######################################################################################
#						    										 				  #
# Spatio-temporal data set with rural background PM10 concentrations in Germany 2005  #
#                           														  #
#######################################################################################

# Set working directory to folder "JCGS_Codes" folder.
setwd("JCGS_Codes")
source("Functions.R")
library(salso)
library(drpm)
library(MCMCpack)
library(mclust)

#
# This data is found in the gstat package
library(gstat)
data(DE_RB_2005)
dat <- DE_RB_2005

# Create ymat with columns corresponding to time rows stations
N <- length(dat@sp)
Tm <- 365
y <- matrix(NA, nrow=N, ncol=Tm)
for(i in 1:Tm){
	y[dat@index[dat@index[,2]==i,1], i] <- dat@data[dat@index[,2]==i,1]
}

# Try to create an average PM10 per month
year <- c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),rep(5,31),
         rep(6,30),rep(7,31),rep(8,31),rep(9,30),rep(10,31),
         rep(11,30),rep(12,31))
week <- rep(1:52, each=7)
ymn <- t(apply(y, 1, function(x) tapply(x,year,mean, na.rm=TRUE)))


## Keep those that don't have any missing values when overageing over a month
ysub2 <- ymn[-c(4,16,25,27,30,43,52,59,69),]


mn <- apply(ysub2,2,mean)
sd <- apply(ysub2,2,sd)

# Center the observations
y <- t(t(ysub2) - mn)


tps <- ncol(y)

s_coords <- (dat@sp@coords)[-c(4,16,25,27,30,43,52,59,69),]
smn <- apply(s_coords,2,mean)
ssd <- apply(s_coords,2,sd)
s_std <- t((t(s_coords) - smn)/ssd)



# modelPriors = c(m0, s20, Asig, Atau, Alam, a, b, be)
# m0 - mean of theta (or phi0)
# s20 - variance of theta (or phi0)
# Asig - maximum value for sigma
# Atau - maximum value for tau
# Alam - maximum value for lambda
# a - shape 1 for alpha
# b - shape 2 for alpha
# b - scale for eta1
modelPriors <- c(0,100, 10, 5, 5, 2, 2, 1)

# m, k0, nu0, L0
cParms <- c(0,1,5,1)

# SIG, TAU, LAM, ETA1, PHI1
mh <- c(1,1, 1, 0.1, 0.1)
sp <- 4

niter=50000; nburn=10000; nthin=40
nout <- (niter-nburn)/nthin

alpha = 0.0
set.seed(1)
models.out <- list()
hh <- 1
# h <- "111"; s <- "0";
model <- "AR1"
for(s in c("0","1")){
	for(h in c("111","110","101","100","011","010","001","000")){

	  m.n <- as.numeric(strsplit(h, "")[[1]])

	  eta1Update <- m.n[1]!=0
	  phi1Update <- m.n[2]!=0
	  alphaUpdate <- m.n[3]!=0

	  if(s=="0"){
		  sc <- NULL
	  } else {
		  sc <- s_std
	  }


	  cat("model is ", h, "\n")
	  cat("space is ", s, "\n")
	  cat("seed is ", 1*hh, "\n")
	  set.seed(1*hh)
	  print(date())
	  out <- drpm_fit(draws=niter, burn=nburn, thin=nthin, y=y, M=1,s_coords=sc,
					          global_alpha=FALSE, modelPriors=modelPriors,
		                alpha_0 = ifelse(alphaUpdate, 0, 1),
		                eta1_0 = ifelse(eta1Update, 0, 1),
		                phi1_0 = ifelse(phi1Update, 0, 1),
					          SpatialCohesion=sp, cParms=cParms,mh=mh)
	  print(date())
    cat("lpml = ", out$lpml, "\n")
    cat("waic = ", out$waic, "\n\n\n")
		models.out[[hh]] <-  out
		names(models.out)[hh] <- paste0("out",h,"_",s,"_",model)

		rho <- list()
		ccprob <- list()

		for(k in 1:tps){
			rho[[k]] <- salso(t(out$Si[k,,]), loss="VI")
		}

		amn <- round(apply(models.out[[hh]]$alpha,2,mean),2)

        # If there is desire to produce plot of each fit uncomment this lines
#		pdf(paste0("PM10_", h,"_",s,"_",model,"_SC",sp,"_2.pdf"),
		    height=10, width=12.5)

			pchs <- c(letters, paste0(letters,0:9))

			par(mfrow=c(3,4))

			for(jj in 1:tps){

				cex1 <- ((y[,jj]-mean(y[,jj]))/sd(y[,jj])+3)/3
				plot(s_std, col=rho[[jj]], pch=pchs[rho[[jj]]],cex=cex1,
							main=bquote(alpha==.(amn[jj]) ~~ Time ~ .(jj)),
							ylab="", xlab="")
			}

#		dev.off()
		hh <- hh + 1
	}
}

# My run of output has been saved in an .RData object.
# It is provided and can be loaded if so desired. 
# load("PM10_ALL_SpCo_4_rev.RData")

# Create table that contains the WAIC and LPML values
lpml <- lpmlr <-  waic <- numeric()
amn <- matrix(NA, nrow = length(models.out), ncol=ncol(y))
for(j in 1:length(models.out)){
	lpml[j] <- models.out[[j]]$lpml
  lpmlr[j] <- lpml.robust(models.out[[j]]$llike)[5]
	waic[j] <- models.out[[j]]$waic
	amn[j,] <- round(apply(models.out[[j]]$alpha,2,mean),2)
}
res <- data.frame(names=names(models.out), lpml=lpml,lpmlr=lpmlr, waic=waic)
res[order(res[,2]),]


library(xtable)
# Order as seen in Table 3 that focuses on Temporal models only
ord <- c(8,7,4,2,7,5,3,1,16,14,12,10,15,13,11,9)
xtable(res[8:1,c(3,4)], digits=0)




library(fields)
# Compute the lagged partitions for all 16 models
ARImats <- vector("list", length(models.out))
for (i in 1:length(models.out)) {
  ARImats[[i]] <- matrix(NA, nrow=tps, ncol=tps)
}

for(h in 1:length(models.out)){
	for(k in 1:tps){
		rho[[k]] <- salso(t(models.out[[h]]$Si[k,,]), loss="binder")
	}

	for(k in 1: tps){
		for(kk in 1: tps){
			ARImats[[h]][k,kk] <- adjustedRandIndex(rho[[k]], rho[[kk]])
		}
	}
}

pch=as.character(1:nrow(y))


#
# This is Figure 5 in the paper.  Includes lagged ARI plots for models that include time only (no space)
#
ord2 <- c(8,6,4,2,7,5,3,1)

# pdf(paste0("LaggedARI_NoSpace2.pdf"), height=7, width=13)

	par(mfrow=c(2,4),mar=c(2,2,2,2), mgp=c(1.5,0.5,0))
	for(h in ord2){
		m.n <- as.numeric(do.call(c, strsplit(strsplit(names(models.out)[h], "t|\\_")[[1]][c(2,3)], "")))
		header <- bquote(eta[1]*.(ifelse(m.n[1]==1,"-Yes","-No"))*","~
	    	             phi[1]*.(ifelse(m.n[2]==1,"-Yes","-No"))*","~
	        	         alpha[t]*.(ifelse(m.n[3]==1,"-Yes","-No")))
		image.plot(ARImats[[h]], main=header, zlim=range(do.call(c,ARImats),na.rm=TRUE), axes=FALSE)
		mtext(text=c(paste("",1:tps)), side=2, line=0.3, at=seq(0,1,length=12), las=1, cex=0.8)
		mtext(text=c(paste("",1:tps)), side=1, line=0.3, at=seq(0,1,length=12), las=2, cex=0.8)
	}
# dev.off()


#  This plot is Figure 6 in the paper.
#  In includes estimated partitions over time for four models.
#
#  temperal dependence in data model 1-Yes, 0-No (eta parameters)
#  temperal dependence in atoms 1-Yes, 0-No (phi1 parameter)
#  temperal dependence in partition 1-Yes, 0-No (alpha parameter)
#
#  spatial dependence in partition
#
# In the plot I include the following models
# 000_0, 001_0, 110_0, 111_0
# pdf("PartitionsOverTime3.pdf", height=15, width=15)
par(mfrow=c(2,2))
for(kkk in c(8,7, 2, 1)){

  rho <- list()
  ccprob <- list()

  for(k in 1:tps){
      rho[[k]] <- salso(t(models.out[[kkk]]$Si[k,,]), loss="binder")
  }

  m.n <- as.numeric(do.call(c, strsplit(strsplit(names(models.out)[kkk], "t|\\_")[[1]][c(2,3)], "")))
  header <- bquote(eta[1]*.(ifelse(m.n[1]==1,"-Yes","-No"))*","~
    	       phi[1]*.(ifelse(m.n[2]==1,"-Yes","-No"))*","~
        	   alpha[t]*.(ifelse(m.n[3]==1,"-Yes","-No")))

  plot(rep(1:12, each=nrow(y)), rep(1:60, times=ncol(y)), type='n',
       yaxt="n",ylab="", xlab="time", main=header, cex.main=2, cex.lab=1.75,
       xaxt="n")
  axis(side=1, at=1:12, cex.axis=1.5)
  for(jj in 1:ncol(y)){
    ord <- order(rho[[jj]])
    cex1 <- rep(0.8,length(ord))
    text(rep(jj, each=nrow(y)), 1:60, labels=pch[ord], col=rho[[jj]][ord], cex=cex1[ord])
    if(jj>1)text(jj-0.5, 61, bquote(alpha[.(jj)]==.(amn[kkk,jj])), cex=0.85)
  }
}

# dev.off()



###########################################
#
# Analysis with space in the partition model
#
#
# This is Figure 7 in the paper.  Lagged-ARI plots for models that include space in the partition
#
ord2 <- c(16,14,12,10,15,13,11,9)

# pdf(paste0("LaggedARI_OnlySpace2.pdf"), height=7, width=13)

par(mfrow=c(2,4),mar=c(2,2,2,2), mgp=c(1.5,0.5,0))
for(h in ord2){
	m.n <- as.numeric(do.call(c, strsplit(strsplit(names(models.out)[h], "t|\\_")[[1]][c(2,3)], "")))
	header <- bquote(eta[1]*.(ifelse(m.n[1]==1,"-Yes","-No"))*","~
    	             phi[1]*.(ifelse(m.n[2]==1,"-Yes","-No"))*","~
        	         alpha[t]*.(ifelse(m.n[3]==1,"-Yes","-No")))
	image.plot(ARImats[[h]], main=header, zlim=range(do.call(c,ARImats),na.rm=TRUE), axes=FALSE)
	mtext(text=c(paste("",1:tps)), side=2, line=0.3, at=seq(0,1,length=12), las=1, cex=0.8)
	mtext(text=c(paste("",1:tps)), side=1, line=0.3, at=seq(0,1,length=12), las=2, cex=0.8)
}
# dev.off()



# This is Figure 8 in the paper.
# It shows spatially refrenced estimated partitions for model 9 (i.e., 111 1 everything turned on)
# which fits data the best.
#
# There are more clusters with sPPM and so I need to create a color vector so that color doesn't
# get replicated and it is easy to identify the clusters

colors <- c("black","red","green","blue","cyan","purple","yellow","gray",
            "orange","maroon", "lavender", "magenta","turquoise")

rho <- list()

for(k in 1:tps){
	rho[[k]] <- salso(t(models.out[[9]]$Si[k,,]), loss="binder")
}
adj.rand <- matrix(NA,nrow=tps, ncol=tps)
for(k in 1:tps){
	for(kk in 1:tps){
		adj.rand[k,kk] <- adjustedRandIndex(rho[[k]], rho[[kk]])
	}
}

# pdf(paste0("PM10_", 111,"_",1,"_",model,"_SC",sp,"4.pdf"), height=12, width=9)

	par(mfrow=c(4,3),mar=c(3,3,2,1), mgp=c(1.5,0.5,0))

	for(jj in 1:tps){
      cex1 <- ((y[,jj]-mean(y[,jj]))/sd(y[,jj])+3)/3
      plot(s_coords, col=rho[[jj]], pch=pchs[rho[[jj]]],cex=cex1,
				main=bquote(Time~period ~ .(jj)*"," ~~ alpha==.(amn[9,jj])),
				ylab="", xlab="", yaxt="n", xaxt="n", type='n')

      text(s_coords[,], labels=pch, col=colors[rho[[jj]]], cex=cex1[])
	}

# dev.off()




#######################################################################
#######################################################################
#
# Table of all adjusted rand index values.  Not included in the paper
#
xtable(round(adj.rand,2))


#
# Plot of adjusted rand index values, along with posterior means of alpha, and overall data values
# This plot was not included in the paper.
#
# pdf("LaggedARIvalues.pdf",height=11, width=9)
par(mfrow=c(4,3))
plot(mn, type='n', xlab="Month", ylab="PM10", main="PM10 for each station",
		ylim=range(t(t(y) + mn)))
for(k in 1:nrow(y)){
	lines(t(t(y) + mn)[k,], type='b' )
}
plot(amn[9,], type='b', xlab="Month", ylab=expression(E(alpha[t]*'|'*Y)), main=expression("Posterior mean of" ~ alpha[t]))
for(k in 1:(tps-2)){
	plot(adj.rand[k,((k+1):tps)], col=1, type='b', xlab="Lagged Month", ylab="Adjusted Rand Index",
#			main= bquote(ARI(hat(rho)[.(k)],.)))
		main= bquote("Lagged ARI values for" ~ hat(rho)[.(k)]))
}
# dev.off()



#
# This plot is not included in the Final draft of the paper
#
# pdf("MeanAlphaPM10.pdf",height=5, width=7)
par(mfrow=c(2,1))
plot(mn, type='b', xlab="Month", ylab="PM10")
plot(amn[9,], type='b', xlab="Month", ylab=expression(E(alpha[t]*'|'*Y)))
# dev.off()


#
# All lagged plots in one figure.  This plot didn't make it into the final draft of the paper.
#
ord <- c(8,6,4,2,7,5,3,1,16,14,12,10,15,13,11,9)

# pdf(paste0("LaggedARI.pdf"), height=10, width=10)

par(mfrow=c(4,4),mar=c(2,2,2,2), mgp=c(1.5,0.5,0))
for(h in ord){
	m.n <- as.numeric(do.call(c, strsplit(strsplit(names(models.out)[h], "t|\\_")[[1]][c(2,3)], "")))
	header <- bquote(eta[1]*.(ifelse(m.n[1]==1,"-Yes","-No"))*","~
    	             phi[1]*.(ifelse(m.n[2]==1,"-Yes","-No"))*","~
        	         alpha[t]*.(ifelse(m.n[3]==1,"-Yes","-No"))*","~
            	     "space"*.(ifelse(m.n[4]==1,"-Yes","-No")))
	image.plot(ARImats[[h]], main=header, zlim=range(do.call(c,ARImats),na.rm=TRUE), axes=FALSE)
	mtext(text=c(paste("",1:tps)), side=2, line=0.3, at=seq(0,1,length=12), las=1, cex=0.8)
	mtext(text=c(paste("",1:tps)), side=1, line=0.3, at=seq(0,1,length=12), las=2, cex=0.8)
}
# dev.off()



