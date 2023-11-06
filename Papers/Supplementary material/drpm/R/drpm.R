#-----------------------------------------------------------------------------------------------------


# R-wrapper that fits the global model with a number of options.
# 1 - temporal dependence in the likelihood
# 2 - temporal dependence in the latent level
# 3 - temporal dependenc in the partition model
# 4 - Spatial information in the partition model.

drpm_fit <- function(y,s_coords=NULL, M=1,
					alpha_0=FALSE,eta1_0=FALSE,phi1_0=FALSE,global_alpha=TRUE,
					modelPriors=c(0,100^2,1,1,1,1,1,1),
					SpatialCohesion=4,
					cParms=c(0, 1, 2, 1),
					mh=c(0.5, 1, 0.1, 0.1, 0.1),
					verbose=FALSE,
					draws=1100,burn=100,thin=1){

	nout <- (draws-burn)/thin

  ntime = ncol(y)
	nsubject = nrow(y)

	s1 <- s_coords[,1]; s2 <- s_coords[,2]

	Si <- llike <- fitted <- matrix(1, nrow=nout, ncol=ntime*nsubject)
	gamma <- matrix(0, nrow=nout, ncol=ntime*nsubject)
	mu <- sig2 <- matrix(1, nrow=nout, ncol=ntime*nsubject)
	theta <- tau2 <- alpha_out <- matrix(0.5, nrow=nout, ncol=ntime);
	eta1 <- matrix(0, nrow=nout, ncol=nsubject);
	phi0 <- phi1 <- lam2 <- rep(0, nout)
	lpml <- waic <- rep(0,1)

	sPPM <- ifelse(is.null(s1), FALSE, TRUE)

	cat("sPPM = ", sPPM, "\n")

#if(!sPPM){
#	C.out <- .C("mcmc_drpm_ar1",
#	C.out <- .C("drpm_ar1_sppm",
#	C.out <- .C("drpm_ar1_sppm2",
#      as.integer(draws), as.integer(burn),as.integer(thin),
#      as.integer(nsubject),as.integer(ntime),
#      as.double(t(y)), as.double(s1), as.double(s2),
#      as.double(M),
#      as.double(modelPriors),as.integer(global_alpha),
#      as.integer(alpha_0), as.integer(eta1_0), as.integer(phi1_0),
#      as.integer(sPPM), as.integer(SpatialCohesion), as.double(cParms),
#      as.double(mh), as.integer(verbose),
#      Si.draws = as.integer(Si),mu.draws = as.double(mu),
#      sig2.draws = as.double(sig2), eta1.draws = as.double(eta1),
#      theta.draws = as.double(theta), tau2.draws = as.double(tau2),
#      phi0.draws = as.double(phi0), phi1.draws = as.double(phi1),
#      lam2.draws = as.double(lam2), gamma.draws=as.integer(gamma),
#      alpha.draws = as.double(alpha_out), fitted.draws = as.double(fitted),
#      llike.draws=as.double(llike), lpml.out = as.double(lpml),
#      waic.out = as.double(waic))
#}


#if(sPPM){
	space_1 <- FALSE;
	alpha <- 0.0
	update_alpha <- ifelse(alpha_0==TRUE, 0, 1)
	update_eta1 <- ifelse(eta1_0==TRUE, 0, 1)
	update_phi1 <- ifelse(phi1_0==TRUE, 0, 1)
	C.out <- .C("drpm_ar1_sppm",
      as.integer(draws), as.integer(burn),as.integer(thin),
      as.integer(nsubject),as.integer(ntime),
      as.double(t(y)), as.double(s1), as.double(s2),
      as.double(M), as.double(alpha),
      as.double(modelPriors),as.integer(global_alpha),
      as.integer(update_alpha), as.integer(update_eta1), as.integer(update_phi1),
      as.integer(sPPM), as.integer(SpatialCohesion), as.double(cParms),
      as.double(mh), as.integer(space_1),
      Si.draws = as.integer(Si),mu.draws = as.double(mu),
      sig2.draws = as.double(sig2), eta1.draws = as.double(eta1),
      theta.draws = as.double(theta), tau2.draws = as.double(tau2),
      phi0.draws = as.double(phi0), phi1.draws = as.double(phi1),
      lam2.draws = as.double(lam2), gamma.draws=as.integer(gamma),
      alpha.draws = as.double(alpha_out),fitted.draws = as.double(fitted),
      llike.draws=as.double(llike),lpml.out = as.double(lpml),
      waic.out = as.double(waic))
#}

	out <- NULL

  out$Si <- array(C.out$Si.draws, c(ntime,nsubject,nout))
  out$gamma <- array(C.out$gamma.draws, c(ntime,nsubject,nout))
  out$mu <- array(C.out$mu.draws, c(ntime,nsubject,nout))
  out$sig2 <- array(C.out$sig2.draws, c(ntime,nsubject,nout))
  out$theta <- matrix(C.out$theta.draws, nrow=nout, byrow=TRUE)
  out$tau2 <- matrix(C.out$tau2.draws,nrow=nout, byrow=TRUE)
  out$alpha <- matrix(C.out$alpha.draws,nrow=nout, byrow=TRUE)
	if(global_alpha) out$alpha <- out$alpha[,2]

  out$eta1 <- matrix(C.out$eta1.draws,nrow=nout, byrow=TRUE)

	out$phi0 <- matrix(C.out$phi0.draws,nrow=nout, byrow=TRUE)
  out$phi1 <- matrix(C.out$phi1.draws,nrow=nout, byrow=TRUE)
  out$lam2 <- matrix(C.out$lam2.draws,nrow=nout, byrow=TRUE)

  out$llike <- matrix(C.out$llike.draws,nrow=nout, byrow=TRUE)
  out$fitted <- matrix(C.out$fitted.draws,nrow=nout, byrow=TRUE)
  out$lpml <- C.out$lpml.out
  out$waic <- C.out$waic.out

  out

}





#-----------------------------------------------------------------------------------------------------



