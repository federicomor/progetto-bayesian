#A <- c(0.25, 0.5, 0.75, 1, 5, 10, 0.5*sd(Yvec))
#waic.mat <- matrix(NA, nrow=5, ncol=length(A))
#for(iii in 1:length(A)){
#    modelPriors <- c(0, 100, A[iii], 100,  100,  1,  1,  1)
#    print(modelPriors)
#    for(iiii in 1:5){
#      cat("iiii = ", iiii, "\n")
#      drpm1 <- drpm_fit(y=Ymat, global_alpha=FALSE,
#                  alpha_0 = FALSE,
#                  eta1_0 = TRUE,
#                  phi1_0 = TRUE,
#                  modelPriors=modelPriors,
#                  draws=niter, burn=nburn, thin=nthin,verbose=TRUE)
#      waic.mat[iiii, iii] <- drpm1$waic
#    }
#}
