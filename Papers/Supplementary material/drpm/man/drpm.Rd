\name{drpm}
\alias{drpm}
\title{Function that fits spatio-temporal hierarchical Gaussian model that is based on a dependent random partition model}
\usage{
drpm(y,
     s_coords=NULL, sp_coords=NULL,
     M=1,
     priorvals,
		 alpha,
		 global_alpha=TRUE,
		 update_alpha=FALSE,
		 update_eta1=FALSE,
		 update_phi1=FALSE,
		 SpatialCohesion=4,
		 cParms=c(0, 1, 2, 1),
		 mh=c(0.5, 1, 0.1, 0.1, 0.1),
		 LikeModel="AR1",
		 draws,
		 burn,
		 thin,
		 verbose=FALSE)
}
\description{
\code{drpm} is the main function used to fit model with Guassian data model and spatio-temporal dependent random partition model.  In the model, temporal dependence could possible appear in three places a) AR(1) process in data model, b) AR(1) process associated with mean model of atmos, c) temporal dependence in the partition.  Details follow regarding how to flip each of these switches on.

}

\arguments{
\item{y}{An m x T numeric matrix containing the response values measured over time.}

\item{s_coords}{Two-column matrix containing spatial locations (i.e., longitude and lattitude).}

\item{sp_coords}{Two-column matrix containing spatial locations at which out-of-sample predictions will be collected.}

\item{M}{Parameter related to Dirichlet process scale or dispersion parameter.}

\item{alpha}{Value assigned to temporal dependence parameter if it is assumed known}

\item{cohesion}{Scalar that indicates which cohesion to use.
\itemize{
    \item 1 - distance from centroids \cr
	  \item 2 - upper bound \cr
	  \item 3 - auxiliary similarity \cr
	  \item 4 - double dipper similarity \cr
  }
}


\item{global_alpha}{Logical indicating if alpha is constant across time.}

\item{update_alpha}{Logical indicating if alpha is assumed zero (i.e., iid partition model)}

\item{update_eta1}{Logical indicating if atoms are modeled with AR(1) or iid}

\item{update_phi1}{Logical indicating if data model is AR(1) or conditionally independent}

\item{SpatialCohesion}{Logical indicating if data model is AR(1) or conditionally independent}


\item{modelPriors}{Vector containing model prior values (see below for more details)}

\item{cParms}{Vector containing partition model prior values (see below for more details)}

\item{mh}{Tuning standard deviations for metropolis updates for sigma2 and sigma20}

\item{draws}{Number of MCMC samples to collect}

\item{burn}{Number of the MCMC samples discarded in the burn-in phase of the sampler}

\item{thin}{The amount of thinning desired for the chain}

}
\details{

The vector \verb{modelPriors = c(m0, s20, ms, ms0)} where each prior parameter is listed in the model description below.

The cParm vector contains values associated with the cohesion function.\cr
\verb{
cParm = c(
  epsilon value - cohesion 1 only,
  distance bound - cohesion 2 only,
  mu0 - center of NNIG for cohesion 3 and 4
  k0 - scale parm of gaussian part of NNIG for cohesion 3 and 4
  v0 - degrees of freedom IG part of NNIG for cohesion 3 and 4
  L0 - scale parm (scalar of identity matrix) IG part of NNIG for cohesion 3 and 4).}

The model this function fits is Gaussian likelihood model using the sPPM prior on partitions (Page and Quintana, 2016).  Specific model details are
  \deqn{y_i | \mu^*, \sigma^{2*}, c_i \sim N(\mu_{c_i}^*, \sigma^{2*}_{c_i}), i=1,\ldots,n}{y_i ~ N(mu*_{c_i}, sigma2*_{c_i}), i=1,\ldots,n}
  \deqn{\mu^*_j | \mu_0, \sigma^2_0 \sim N(\mu_0,\sigma^2_0)}{mu*_j | mu0, sig20 ~ N(mu0,sig20)}
  \deqn{\sigma^*_j | A \sim UN(0,  ms)}{sigma*_j | A  ~ UN(0, ms)}
  \deqn{\rho|M,\xi  \sim sPPM}{rho ~ sPPM}

To complete the model specification, the folloing hyperpriors are assumed,
  \deqn{\mu_0 | m, s^2 \sim N(m0,s0^2)}{mu0 | m0, s02 ~ N(m0,s02)}
  \deqn{\sigma_0 | B  \sim UN(0,ms0)}{sigma0 | B ~ UN(0,ms0)}

Note that we employ uniform prior distributions on variance components as suggest in Gelman's 2006 Bayesian paper.  "sPPM" in the model specificaiton denotes the the spatial product partition model.  The computational implementation of the model is based algorithm 8 found in Neal's 2000 JCGS paper.
}

\value{
This function returns in a list all MCMC interates for each model parameter, posterior predictive, and fitted values.  In addition the LPML model fit metric is provided.
}
\examples{


\dontrun{

# library(akima) #interp function
# library(maps) #maps function
# library(geoR)
# library(salso) # to estimate partition
# library(ppmSuite)

data(scallops)

Y<-log(scallops[,5]+1)
s_coords <- scallops[,3:4] #lat and long
m <- dim(s_coords)[1]


# As an exloratory data analysis plot data
# and fit a variogram

# plot data and fit contour lines using interpolator
map("usa", xlim=c(-75.5, -71), ylim=c(38, 42))
cex1 <- ((Y-mean(Y))/sd(Y)+3)/3
points(s_coords[,2], s_coords[,1], cex=cex1, pch=1)

int.scp <- interp(s_coords[,2], s_coords[,1], Y) ## interpolator need akima package
contour(int.scp, add=TRUE)

# Variogram on raw data
# Arrange the data in geoR format:
geodat <- list(coords=cbind(s_coords[,2],s_coords[,1]),data=Y)

# Compute the empirical variogram using variog
bins<-seq(0,2.5,length=20)  #there are the bin edges
vg <- variog(geodat, uvec=bins)
plot(vg)

fit2 <- variofit(vg,cov.model="exponential",ini.cov.pars=c(3,1))
lines(fit2, lty=3, col='blue')

# standardize spatial coordinates
smn <- apply(s_coords,2,mean)
ssd <- apply(s_coords,2,sd)
s_std <- t((t(s_coords) - smn)/ssd)

# Create a grid of prediction locations
np <- 100
sps <- interp(x=s_coords[,1], y=s_coords[,2], z=Y, xo=seq(min(s_coords[,1]), max(s_coords[,1]), length = np),
       yo=seq(min(s_coords[,2]), max(s_coords[,2]), length = np))

sp1 <- expand.grid(seq(min(s_coords[,1]), max(s_coords[,1]),length=np), seq(min(s_coords[,2]), max(s_coords[,2]), length=np))
## Keep only points in convex hull
sp <- sp1[!is.na(c(sps$z)),]

sp_std <- t((t(sp) - smn)/ssd) # standardized prediction spatial coordinates


niter <- 20000
nburn <- 10000
nthin <- 10
nout <- (niter - nburn)/nthin


out <- sppm(y=Y,s=s_std,s.pred=sp_std,Cohesion=4, M=1, draws=niter, burn=nburn, thin=nthin)

# fitted values
fitted.values <- out$fitted
fv.mn <- apply(fitted.values, 2,mean)
mean((Y - fv.mn)^2) # MSE
out$lpml #lpml value

ppred <- out$ppred
predmn <- apply(ppred,2,mean)

qplot(sp[,1], sp[,2], geom="tile", fill=predmn)

# estimate and plot the partition
nclus <- out$nclus
plot(nclus, type='l')
Si <- out$Si
probs <- psm(Si, parallel=FALSE)
part.est <- salso(probs, loss="binder", nPermutations=50, parallel=FALSE)$estimate
nc <- length(unique(part.est)); nk <- table(part.est)
plot(s_coords, cex=cex1,pch=letters[part.est], col=tim.colors(nc)[part.est], ylab="latitude", xlab="longitude")

# plot the number of clusters
nclus <- out$nclus
mean(nclus)
plot(nclus, type='l', xlab="MCMC iterate", ylab="k")
plot(table(nclus)/nout, type='h', ylab="Pr(k)", xlab="Number of Clusters")




  }


}
