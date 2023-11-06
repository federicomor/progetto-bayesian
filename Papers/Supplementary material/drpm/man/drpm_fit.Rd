\name{drpm_fit}
\alias{drpm_fit}
\title{Function that fits spatio-temporal hierarchical Gaussian model that is based on a dependent random partition model}
\usage{
drpm_fit(y,
         s_coords=NULL,
         M=1,
         alpha_0=FALSE,
         eta1_0=FALSE,
         phi1_0=FALSE,
         global_alpha=TRUE,
         modelPriors=c(0,100^2,1,1,1,1,1,1),
         SpatialCohesion=4,
         cParms=c(0, 1, 2, 1),
         mh=c(0.5, 1, 0.1, 0.1, 0.1),
         verbose=FALSE,
         draws=1100,
         burn=100,
         thin=1)
}


\description{
\code{drpm_fit} is the main function used to fit model with Guassian data model and spatio-temporal dependent random partition model.  In the model, temporal dependence could possible appear in three places a) AR(1) process in data model, b) AR(1) process associated with mean model of atmos, c) temporal dependence in the partition.  Details follow regarding how to flip each of these switches on.

}

\arguments{
\item{y}{An m x T numeric matrix containing the response values measured over time.}

\item{s_coords}{Two-column matrix containing spatial locations (i.e., longitude and lattitude).}

\item{M}{Parameter related to Dirichlet process scale or dispersion parameter.}

\item{global_alpha}{Logical indicating if alpha is constant across time.}

\item{alpha_0}{Logical where TRUE indicates alpha_t = 0 (i.e., iid partition model used)}

\item{eta1_0}{Logical where TRUE indicates eta1 = 0 (i.e., observations are conditionally independent).}

\item{phi1_0}{Logical where TRUE indicates phi1 = 0 (i.e., iid model for atoms)}

\item{modelPriors}{Vector containing the following hierarchical model prior values in order
\itemize{
  \item {m0 - mean of phi0}
	\item {s20 - variance of phi0}
	\item {A - upper bound on sigma }
	\item {At - upper bound on tau }
	\item {Al - Upper bound on lambda }
	\item {at - shape 1 beta prior parameter for alpha }
	\item {bt - shape 2 beta prior parameter for alpha }
	\item {be - scale prior for xi}
  }
}
\item{SpatialCohesion}{Logical indicating what cohesion to use if spatial coordinates are provided
\itemize{
\item {3 - auxiliary similarity }
\item {4 - double dipper similarity }
  }
}


\item{cParms}{Vector containing cohesion parameter values (see below for more details)
\itemize{
\item {mu0 - center of NNIG }
\item {k0 - scale parameter of gaussian part of NNIG }
\item {v0 - degrees of freedom for inverse-gamma part of NNIG}
\item {L0 - scale parameter for inverse-gamma part of NNIG}
  }
}
\item{mh}{Vector of gaussian standard deviations for metropolis updates of sigma2, tau, lambda, eta1, phi1.}

\item{draws}{Number of MCMC samples to collect}

\item{burn}{Number of the MCMC samples discarded in the burn-in phase of the sampler}

\item{thin}{The amount of thinning desired for the chain}

\item{verbose}{Logical indicating if MCMC progression should be print to screen}
}

\value{
This function returns, in a list, all MCMC interates for each model parameter.  In addition the LPML and WAIC model fit metric is provided.
}
\examples{


\dontrun{


   Still needs to be done.

  }


}
