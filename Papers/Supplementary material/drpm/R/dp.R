gaussian_dp <- function(y, N, m, v, A, A0, mh=c(0.5, 0.5), alpha, niter, nburn, nthin){

  n <- length(y)

  run <- .Call("GIBBS",
          as.double(y), as.integer(n), as.integer(N),
          as.double(m), as.double(v), as.double(A), as.double(A0),
          as.double(alpha), as.double(mh),
          as.integer(niter), as.integer(nburn),
				  as.integer(nthin))

  run
}
