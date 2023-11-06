## I am standardizing the x's.  For
linearDDP <- function(draws,burn,thin,y,x,ygrid,xpred,priorparms,Mdp){

	# prior

	out <- NULL

	nobs <- length(y)
	nygrid <- length(ygrid)
	nxpred <- ifelse(is.null(xpred), 1, length(xpred))

	nout <- (draws - burn)/thin

	Mval <- Mdp
	cat("M = ", Mval, "\n")
	X <- cbind(rep(1,nrow(x)),x)
	ncov <- ncol(X)

	cat("ncov = ", ncov, "\n")

	Xpred <- cbind(rep(1, nxpred),xpred)

	density <- rep(0,nout*nygrid*nxpred)
	cumdensity <- rep(0,nout*nygrid*nxpred)
	beta <- rep(0,nout*nobs*ncov)
	sig2 <- rep(0,nout*nobs)
	beta0 <- rep(0,nout*ncov)
	Sig0 <- rep(0,nout*ncov*ncov)
	Si <- rep(1,nout*nobs)
	nclus <- rep(1,nout)
	M <- rep(1, nout); if(Mdp>0) M <- rep(Mdp,1)
	lpml <- rep(0,1)
	waic <- rep(0,1)
	llike <- rep(0,nout*nobs)
	fitted <- rep(0,nout*nobs)

	cat("length Sig0 = ", length(Sig0), "\n")


	C.out <- .C("LDDP",
              	as.integer(draws), as.integer(burn),as.integer(thin),
              	as.integer(nobs),as.integer(ncov),as.double(y), as.double(t(X)),
              	as.double(Mval),as.integer(nygrid), as.double(ygrid),
              	as.integer(nxpred), as.double(t(Xpred)),
              	as.double(priorparms),
              	betah.draws = as.double(beta),
              	sig2h.draws = as.double(sig2),
              	beta0.draws = as.double(beta0),
              	Sig0.draws = as.double(Sig0),
              	M.draws = as.double(M),
              	Si.draws = as.integer(Si),
              	nclus.draws = as.integer(nclus),
              	density.draws = as.double(density),
              	cumdensity.draws= as.double(cumdensity),
              	fitted.draws = as.double(fitted),
              	llike.draws = as.double(llike),
              	lpml.out = as.double(lpml),waic.out = as.double(waic)
              	)


	out$betah <- array(C.out$betah.draws,c(ncol(X), nobs, nout))
	out$sig2h <- matrix(C.out$sig2h.draws,nrow=nout, byrow=TRUE)
	out$beta0 <- matrix(C.out$beta0.draws,nrow=nout, byrow=TRUE)
	out$Sig0 <- matrix(C.out$Sig0.draws,nrow=nout, byrow=TRUE)
	out$M <- matrix(C.out$M.draws,nrow=nout, byrow=TRUE)
	out$Si <- matrix(C.out$Si.draws,nrow=nout, byrow=TRUE)
	out$nclus <- C.out$nclus.draws
 	out$density <- matrix(C.out$density.draws,nrow=nout, byrow=TRUE)
 	out$cumdensity <- matrix(C.out$cumdensity.draws,nrow=nout, byrow=TRUE)
 	out$fitted.values <- matrix(C.out$fitted.draws,nrow=nout, byrow=TRUE)
 	out$llike <- matrix(C.out$llike.draws,nrow=nout, byrow=TRUE)
 	out$lpml <- C.out$lpml.out
 	out$waic <- C.out$waic.out
 	out

}

