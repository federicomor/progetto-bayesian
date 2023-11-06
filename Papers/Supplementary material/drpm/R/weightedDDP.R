# R-wrapper that fits a joint model between y and x and
# employes the conditional to assess relationship between them.
# x can be a covariate or time.

## I am standardizing the x's.  For
weightedDDP <- function(draws,burn,thin,y,x,ygrid,xpred,priorparms,N, Mdp){

	# prior

	out <- NULL

	nobs <- length(y)
	nygrid <- length(ygrid)
	nxpred <- ifelse(is.null(xpred), 1, length(xpred))

	nout <- (draws - burn)/thin

	Mval <- Mdp
	cat("M = ", Mval, "\n")

	cat("N = ", N, "\n")

	Xpred <- cbind(rep(1, nxpred),xpred)


	Ymat <- cbind(y, x)
  dim <- ncol(Ymat)

	Si <- rep(1,nout*nobs)
  mu <- rep(0,nout*N*dim)
	Sigma <- rep(0,nout*dim*dim*N)
	Si <- rep(1,nout*nobs)
	mu0 <- rep(0, nout*dim)
	Sigma0 <- rep(0, nout*dim*dim)
	sbw <- rep(1,nout*N)
  nclus <- rep(1, nout)

	nclus <- rep(1,nout)
	M <- rep(1, nout); if(Mdp>0) M <- rep(Mdp,1)
	lpml <- rep(0,1)
	waic <- rep(0,1)
	llike <- rep(0,nout*nobs)
	fitted <- rep(0,nout*nobs)

	C.out <- .C("WDDP",
              	as.integer(draws), as.integer(burn),as.integer(thin),
              	as.integer(nobs),as.integer(dim), as.integer(N),
	              as.double(t(Ymat)), as.double(priorparms),
              	mu.draws = as.double(mu),
              	Sigma.draws = as.double(Sigma),
              	Si.draws = as.integer(Si),
              	nclus.draws = as.integer(nclus),
	              mu0.draws = as.double(mu0),
	              Sigma0.draws = as.double(Sigma0),
	              sbw.draws = as.double(sbw),
	              llike.draws = as.double(llike),
	              fitted.draws = as.double(fitted),
              	lpml.out = as.double(lpml),
	              waic.out = as.double(waic)
              )


	out$mu <- array(C.out$mu.draws,c(dim, N, nout))
	out$Sigma <- matrix(C.out$Sigma.draws,nrow=nout, byrow=TRUE)
	out$Si <- matrix(C.out$Si.draws,nrow=nout, byrow=TRUE)
	out$sbw <- matrix(C.out$sbw.draws,nrow=nout, byrow=TRUE)
	out$mu0 <- matrix(C.out$mu0.draws,nrow=nout, byrow=TRUE)
	out$Sigma0 <- matrix(C.out$Sigma0.draws,nrow=nout, byrow=TRUE)
	out$llike <- matrix(C.out$llike.draws,nrow=nout, byrow=TRUE)
	out$fitted <- matrix(C.out$fitted.draws,nrow=nout, byrow=TRUE)
	out$nclus <- C.out$nclus.draws
 	out$lpml <- C.out$lpml.out
 	out$waic <- C.out$waic.out
 	out

}


