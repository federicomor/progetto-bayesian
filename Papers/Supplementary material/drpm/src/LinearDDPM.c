/*************************************************************
 * Copyright (c) 2016 Garritt Leland Page
 *
 * This file contains C code for an MCMC algorithm constructed
 * to fit a density using a Dependent Dirichlet Process Mixture (DDPM).
 * I will employ algorithm 8 of Neal with m=1
 *
 * The following density model 
 *
 * f(Y|x) ~ int phi(y|theta) dG_x(theta) 
 * 
 * can be written hierarchically after introducing cluster
 * labels 
 *
 * Y_i | x_i, c_i mu^*, sigma^2* ~ N(mu^*_{c_i}(x_i), sigma^2*_{c_i})
 
 * with Priors:
 *
 * mu_j(x_i),sigma2_j ~ N(mu0, sig20)xIG(asig,bsig)
 * mu0, sig20 ~ N(m,s2) x IG(a0,b0)
 * 
 *************************************************************/


#include "matrix.h"
#include "Rutil.h"

#include <R_ext/Lapack.h>
#include <R.h>
#include <Rmath.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*****************************************************************************************
* The following are the inputs of the function that are read from R
*
* draws = total number of MCMC draws
* burn = number of MCMC draws discarded as burn-in
* thin = indicates how much MCMC chain should be thinned

* nobs = scalar indicating the number of observations in diseased group
* ncov = scalar indicating the number of covariates
* y = nobs x 1 vector that contains response values
* X = nobs x ncov matrix containing the covariates
* 
* Mval = scalar that contains value of DP precision parameter (Mval=0 implies unknown).
* 
* nygrid = integer scalar indicating the number of grid values
* ygrid = vector of ngrid values at which to evaluate the density for latent variable.
* 
* nxpred = integer scalar indicating the number of x values to make predicions 
* Xpred = nxgrid x p containing covariate values at which to predict.
*
* co = vector containing fixed cut-off points for associated with latent variable
*
* priorparms = vector of prior parameter values the entires correspond to he following
* 	0 - prior for the mean of beta0 (q + nknots + 1) dimension with same value
*   1 - prior for the variance of beta0 assuming that variance for each coefficient same
* 	2 - prior on the diagonal of S0 (initial guess of Sig0)
*   3 - prior on the degrees of freedom for the wishart
*   4 - prior shape of the inverse gamma for sig2
*   5 - prior rate of the inverse gamme for sig 2
*
* The remaining arguments are empty vectors that will pass the MCMC iterates back to R
* 
*****************************************************************************************/


void LDDP(int *draws, int *burn, int *thin, int *nobs, int *ncov,
			  double *y, double *X, double *Mval, int *nygrid, double *ygrid, int *nxpred, 
			  double *Xpred, double *priorparms, 
			  double *beta, double *sig2, double *beta0, 
			  double *Sig0, double *M, int *Si, int *nclus,
			  double *density, double *cumdensity, double *fitted,
			  double *llike, double *lpml, double *waic){

	// i - MCMC iterate
	// ii - MCMC iterate that is saved
	// j - observation iterate
	// jj - second observation iterate
	// k - cluster iterate
	// b - b-spline coefficient iterate
	// bb - second b-spline coefficient iterate
	// g - xgrid iterate;
	
	int i, ii, j, jj, k, b, bb;

	Rprintf("nobs = %d\n", *nobs);

	int nout = (*draws-*burn)/(*thin);

	Rprintf("nout = %d\n", nout);

//	RprintVecAsMat("Xcon", Xcon, N, *ncon);
//	RprintIVecAsMat("Xcat", Xcat, N, *ncat);
	Rprintf("Mval = %f\n", *Mval);

	int nb = *ncov;
	Rprintf("nb = %d\n", nb);
//	RprintVecAsMat("y", y, 1, *nobs);
//	RprintVecAsMat("X", X, *nobs, nb);
	

	Rprintf("nygrid = %d\n", *nygrid);
	Rprintf("nxpred = %d\n", *nxpred);

//	RprintVecAsMat("ygrid", ygrid, 1, *nygrid);
//	RprintVecAsMat("Xpred", Xpred, *nxgrid, nb);
	
	// ===================================================================================
	// 
	// Memory vectors to hold MCMC iterates for non cluster specific parameters
	// 
	// ===================================================================================
		
	double *beta0_iter = R_VectorInit(nb, 0.0);
	double *Sig0_iter = R_VectorInit(nb*nb, 0.0);;
	double *Sig0I = R_VectorInit(nb*nb, 0.0);;
	identity_matrix(Sig0_iter, nb);
	identity_matrix(Sig0I, nb);

	double *Sig0Ibeta0 = R_VectorInit(nb,0.0);
	matrix_product(Sig0I, beta0_iter, Sig0Ibeta0, nb, 1, nb);
	
	double M_iter = 1.0;

	if(*Mval > 0) M_iter = *Mval;

	int Si_iter[*nobs];
	int nclus_iter = 0;

	// ===================================================================================
	// 
	// Memory vectors to hold MCMC iterates for cluster specific parameters
	// 
	// ===================================================================================

	double *sig2h = R_VectorInit(*nobs, 0.1);
	double *betah = R_VectorInit((*nobs)*(nb), 0.0);

	

	int nh[*nobs];


	// ===================================================================================
	// 
	// Initialize a cluster labels and cluster size
	// 
	// ===================================================================================

	// Initialize Si with all individuals in one cluster
	for(j = 0; j < *nobs; j++){ 
		Si_iter[j] = 1; 
		nh[j] = 0;
	}


	// Initial enumeration of number of observations per cluster;
	for(j = 0; j < *nobs; j++){
	
		for(k = 0; k < *nobs; k++){
		
			if(Si_iter[j] == k+1) nh[k] = nh[k] + 1;
		}
	}
	// Initialize the number of clusters	
	for(j = 0; j < *nobs; j++){
	
		if(nh[j] > 0) nclus_iter = nclus_iter + 1;
	
	}


//	Rprintf("nclus_iter = %d\n", nclus_iter);
//	RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);
//	RprintIVecAsMat("nh", nh, 1, nclus_iter);

	
	

	// ===================================================================================		
	//
	// scratch vectors of memory needed to update parameters
	//
	// ===================================================================================

	// Stuff I need to update zi
	double s2tmp;

	// stuff I need to update Si
	int iaux;

	double auxs, auxb, uu;
	double sumsq, denph, cprobh, maxph, xb;
	double sigdraw; 
	double *betadraw = R_VectorInit(nb, 0.0);
	
	double *ph = R_VectorInit(*nobs, 0.0);
	double *probh = R_VectorInit(*nobs, 0.0);


	// stuff I need to update sig2h (player specific), beta0, Sig20;
	double astar, bstar, ld, lln, llo, llr;

	
	// stuff I need to update betah, beta0
	double *scr1 = R_VectorInit(nb*nb, 0.0);
	double *scr2 = R_VectorInit(nb*nb, 0.0);
	double *scr3 = R_VectorInit(nb*nb, 0.0);
	double *sumbetah = R_VectorInit(nb, 0.0);
	double *sumyxbeta = R_VectorInit(nb, 0.0);
	double *Mstar = R_VectorInit(nb, 0.0);
	double *Sstar = R_VectorInit(nb*nb, 0.0);
	double *outrmvnorm = R_VectorInit(nb*nb, 0.0);

	// stuff I need to update Sig0
	double nustar;
	double *bstar_b = R_VectorInit(nb, 0.0);
	double *S0star = R_VectorInit(nb*nb, 0.0);
	double *outrwish = R_VectorInit(nb*nb, 0.0);
	double *zerov = R_VectorInit(nb*nb, 0.0);

	// Stuff to compute lpml, likelihood, WAIC, and Rao-Blackwellized density values
	double lpml_iter, elppdWAIC;
	double *CPOinv = R_VectorInit(*nobs, 0.0);
	double *like_iter = R_VectorInit(*nobs, 0.0);
	double *mnlike = R_VectorInit(*nobs, 0.0);
	double *mnllike = R_VectorInit(*nobs, 0.0);



	// ===================================================================================		
	//
	// Prior parameter values
	//
	// ===================================================================================

//	RprintVecAsMat("priorparms", priorparms, 1, 5);


	// priors for beta 
	double *bvec= R_VectorInit(nb,priorparms[0]); 
//	bvec[0] = 0.5;
	double s2=priorparms[1]; 

//	RprintVecAsMat("bvec",bvec,1,nb);
//	Rprintf("s2 = %f\n", s2);


	// priors for Sig0 
	double *S0=R_VectorInit(nb*nb, 0);

	identity_matrix(S0, nb);
	for(b = 0; b < nb*nb; b++) S0[b] = S0[b]*priorparms[2];

//	RprintVecAsMat("S0", S0, nb, nb);


	double nu0=priorparms[3];
//	Rprintf("nu0 = %f\n", nu0);

	// prior values for sigh
	double A=priorparms[4];
	double A0 = 0.0;
//	Rprintf("A = %f\n", A);

	double asig=1.0, bsig=1.0;
	
	
    Rprintf("Prior values: mb0 = %.2f, s2b0 = %.2f\n \t asig = %.2f, bsig = %.2f\n \t S0 = %.2f, nu0 = %.2f\n",
      	           priorparms[0], s2, asig, bsig, priorparms[2], nu0);
	

	GetRNGstate();


	// ===================================================================================
	//
	// start of the mcmc algorithm;
	//
	// ===================================================================================
	ii = 0;
	for(i = 0; i < *draws; i++){

		if((i+1) % 10000 == 0){
			time_t now;
			time(&now);

			Rprintf("mcmc iter = %d =========================================== \n", i+1);
			Rprintf("%s", ctime(&now));
//			RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);
//			RprintVecAsMat("tau2h", tau2h, 1, nclus_iter);
		}


		//////////////////////////////////////////////////////////////////////////////////
		//
		// update the cluster labels using the polya urn scheme of
		// algorithm 8 found in  Radford Neal's 
		//	"Markov Chain Sampling Methods for Dirichlet Process Mixture Models" 
		//	paper.	
		//
		//////////////////////////////////////////////////////////////////////////////////
		
//		RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);
//		Rprintf("nclus_iter = %d\n", nclus_iter);
		for(j = 0; j < *nobs; j++){

		
//			Rprintf("j = %d =================== \n", j);


//			RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);
//			RprintIVecAsMat("nh", nh, 1, *nobs);
//			Rprintf("nclus_iter = %d\n", nclus_iter);

			
			if(nh[Si_iter[j]-1] > 1){
				
				// Observation belongs to a non-singleton ...
				nh[Si_iter[j]-1] = nh[Si_iter[j]-1] - 1;
				
			}else{
			
				// Observation is a member of a singleton cluster ...
				
				iaux = Si_iter[j];
//				Rprintf("iaux = %d\n", iaux);
				if(iaux < nclus_iter){
				
					// Need to relabel clusters. I will do this by swapping cluster labels
					// Si_iter[j] and nclus_iter along with cluster specific parameters;
				
				
					// All members of last cluster will be assigned subject j's label
					for(jj = 0; jj < *nobs; jj++){
					
						if(Si_iter[jj] == nclus_iter){
						
							Si_iter[jj] = iaux; 
						
						}
					
					}
					
				
					Si_iter[j] = nclus_iter;

					// The following steps swaps order of cluster specific parameters
					// so that the newly labeled subjects from previous step retain
					// their correct cluster specific parameters
					auxs = sig2h[iaux-1];
					sig2h[iaux-1] = sig2h[nclus_iter-1];
					sig2h[nclus_iter-1] = auxs;

					for(b = 0; b < nb; b++){

						auxb = betah[(iaux-1)*(nb) + b];
						betah[(iaux-1)*(nb) + b] = betah[(nclus_iter-1)*(nb) + b];
						betah[(nclus_iter-1)*(nb) + b] = auxb;

					}


					// the number of members in cluster is also swapped with the last
					nh[iaux-1] = nh[nclus_iter-1];
					nh[nclus_iter-1] = 1;
			
				}
			
			
				// Now remove the ith obs and last cluster;
				nh[nclus_iter-1] = nh[nclus_iter-1] - 1;
				nclus_iter = nclus_iter - 1;
			
			
			}
			
//			RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);

//			RprintVecAsMat("sig2h", tau2h, 1, *nobs);
//			RprintVecAsMat("muh", lamh, 1, *nobs);
				
//			Rprintf("nclus_iter = %d\n", nclus_iter);

		

			// The atoms have been relabeled if necessary and now we need to 
			// update Si.


			///////////////////////////////////
			//	
			// Begin the cluster probabilities
			//
			//////////////////////////////////
			for(k = 0 ; k < nclus_iter; k++){

//				Rprintf("k = %d ==================== \n", k);
				xb = 0.0;
				for(b = 0; b < nb; b++){

					xb = xb + betah[k*(nb)+b]*X[j*(nb)+b];

				}
				
				
				ph[k] = dnorm(y[j], xb, sqrt(sig2h[k]), 1) +
				        log((double) nh[k]) - 
				        log(*nobs - 1 + M_iter);

//				Rprintf("ph[k] = %f\n", ph[k]);
			}




//			RprintVecAsMat("ph", ph, 1, nclus_iter);

			for(b = 0; b < nb*nb; b++) scr2[b] = Sig0_iter[b];
			cholesky(scr2, nb , &ld);
			ran_mvnorm(beta0_iter, scr2, nb, scr1, betadraw);

//			RprintVecAsMat("beta0_iter", beta0_iter, 1, nb);
//			RprintVecAsMat("Sig0_iter", Sig0_iter, nb, nb);
//			RprintVecAsMat("betadraw", betadraw, 1, nb);

			sigdraw = sqrt(1.0/rgamma(asig,bsig)); // shape and scale  E(sig2) = asigbsig
//			Rprintf("sig2draw = %f\n", sig2draw);

//			sigdraw = runif(A0,A); //

			xb = 0.0;
			for(b = 0; b < nb; b++){

				xb = xb + betadraw[b]*X[j*(nb)+b];

			}


			ph[nclus_iter] = dnorm(y[j], xb, sigdraw, 1) +
			                 log(M_iter) -  
			                 log(*nobs - 1 + M_iter);


//			RprintVecAsMat("ph", ph, 1, nclus_iter+1);

					
			maxph = ph[0];
			for(k = 1; k < nclus_iter+1; k++){
				if(ph[k] > maxph) maxph = ph[k];
			}		
//			Rprintf("maxph = %f\n", maxph);
			
			denph = 0.0;
			for(k = 0; k < nclus_iter+1; k++){
						
				ph[k] = exp(ph[k] - maxph);
//				ph[k] = pow(exp(ph[k] - maxph), (1 - exp(-0.0001*(i+1))));
				denph = denph + ph[k];
					
			}

//			RprintVecAsMat("ph", ph, 1, nclus_iter+1);

			for(k = 0; k < nclus_iter+1; k++){

				probh[k] = ph[k]/denph;

			}
//			Rprintf("denph = %f\n", denph);

//			RprintVecAsMat("probh", probh, 1, nclus_iter+1);
					
			uu = runif(0.0,1.0);
//			Rprintf("uu = %f\n", uu);


		
			cprobh= 0.0;;
			for(k = 0; k < nclus_iter+1; k++){

				cprobh = cprobh + probh[k];

				if (uu < cprobh){
									
					iaux = k+1;
					break;	
				}
			}		


// 			Rprintf("iaux = %d\n \n \n", iaux);

			if(iaux <= nclus_iter){

				Si_iter[j] = iaux;
				nh[Si_iter[j]-1] = nh[Si_iter[j]-1] + 1;

			}else{
			
				nclus_iter = nclus_iter + 1;
				Si_iter[j] = nclus_iter;
				nh[Si_iter[j]-1] = 1;			
				
				for(b = 0; b < nb; b++){
					betah[(Si_iter[j]-1)*nb + b] = betadraw[b];
				}
	
				sig2h[Si_iter[j]-1] = sigdraw*sigdraw;

			}


		}

	
//		RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);
//		Rprintf("nclus_iter = %d\n", nclus_iter);
//		RprintIVecAsMat("nh", nh, 1, nclus_iter);
		




	
//		RprintVecAsMat("sig2h", sig2h,1,nclus_iter);
//		RprintVecAsMat("betah", betah,nclus_iter,nb);



		//////////////////////////////////////////////////////////////////////////////////
		//																				//
		// udpate betah cluster specific means;											//
		//																				//
		//////////////////////////////////////////////////////////////////////////////////


		for(k = 0; k < nclus_iter; k++){
			
//			Rprintf("k = %d ====== \n", k);
			
			for(b = 0; b < nb; b++){
				sumyxbeta[b] = Sig0Ibeta0[b];
				for(bb = 0; bb < nb; bb++){
					Sstar[b*nb+bb] = Sig0I[b*nb+bb];
				}
			}

//			RprintVecAsMat("sumyxbeta", sumyxbeta, 1, nb);
//			RprintVecAsMat("Sstar",Sstar, nb, nb);


			for(j = 0; j < *nobs; j++){
				if(Si_iter[j] == k+1){
			
					for(b = 0; b < nb; b++){
						sumyxbeta[b] = sumyxbeta[b] + (1/sig2h[k])*y[j]*X[j*(nb) + b];
						for(bb = 0; bb < nb; bb++){

							Sstar[b*nb+bb] = Sstar[b*nb+bb] + (1/sig2h[k])*
							 								   X[j*(nb) + b]*
							                                   X[j*(nb) + bb];
						}
					}
				}
			}
	

//			RprintVecAsMat("sumyxbeta", sumyxbeta, 1, nb);
//			RprintVecAsMat("Sstar",Sstar, nb, nb);

			cholesky(Sstar, nb, &ld);
			inverse_from_cholesky(Sstar, scr1, scr2, nb); 

//			RprintVecAsMat("Sstar",Sstar, nb, nb);
//				
			matrix_product(Sstar, sumyxbeta, Mstar, nb, 1, nb);
			
//			RprintVecAsMat("Mstar", Mstar, 1, nb);
//			RprintVecAsMat("Sstar", Sstar, nb, nb);


			cholesky(Sstar, nb , &ld);
															
			ran_mvnorm(Mstar, Sstar, nb, scr1, outrmvnorm);


//			RprintVecAsMat("betah",outrmvnorm, 1, nb);

			for(b = 0; b < nb; b++){
			
				betah[k*nb + b] = outrmvnorm[b]; 
			
			}
			
			
		}

//		RprintVecAsMat("betah", betah, nclus_iter, nb);


	
		//////////////////////////////////////////////////////////////////////////////////
		//																				//
		// udpate beta0 mean of betah     												//
		//																				//
		//////////////////////////////////////////////////////////////////////////////////
		for(b = 0; b < nb; b++) sumbetah[b] = 0.0;

		for(k = 0; k < nclus_iter; k++){
			for(b = 0; b < nb; b++){
			
				sumbetah[b] = sumbetah[b] + betah[k*nb + b];
			}			
		}


		matrix_product(Sig0I,sumbetah, scr1, nb, 1, nb);

		for(b=0; b < nb; b++){
			scr1[b] = scr1[b] + (1/s2)*bvec[b];
			for(bb = 0; bb < nb; bb++){

				Sstar[b*nb+bb] = Sig0I[b*nb+bb]*nclus_iter;

				if(bb==b) Sstar[b*nb+bb] = Sig0I[b*nb+bb]*nclus_iter + (1/s2);	
			}
		}

		cholesky(Sstar, nb, &ld);
		inverse_from_cholesky(Sstar, scr3, scr2, nb); 

		matrix_product(Sstar, scr1, Mstar, nb, 1, nb);

		cholesky(Sstar, nb , &ld);


//		RprintVecAsMat("Mstar", Mstar, 1, nb);
//		RprintVecAsMat("Sstar", Sstar, nb, nb);
															
		ran_mvnorm(Mstar, Sstar, nb, scr1, outrmvnorm);
		
		for(b = 0; b < nb; b++){
			
			beta0_iter[b] = outrmvnorm[b]; 
			
		}
		
//		RprintVecAsMat("beta0_iter", beta0_iter, 1, nb);



		//////////////////////////////////////////////////////////////////////////////////
		//																				//
		// udpate Sig0 variance matrix associated with betah     			            //
		//																				//
		//////////////////////////////////////////////////////////////////////////////////
		for(b=0; b < nb*nb; b++) S0star[b] = S0[b];
	
	
		for(k = 0; k < nclus_iter; k++){
			for(b = 0; b < nb; b++){

				bstar_b[b] = (betah[k*nb+b] - beta0_iter[b]);
			
			}
	

			matrix_product(bstar_b, bstar_b, scr1, nb, nb, 1);



			for(bb = 0; bb < nb*nb; bb++){
				S0star[bb] = S0star[bb] + scr1[bb];		
			}
						


		}

					
		cholesky(S0star, nb, &ld);
					
		// I need to invert Astar to get an inverse wishart from the ran_wish function;
		inverse_from_cholesky(S0star, scr1, scr2, nb); 
					
		cholesky(S0star, nb, &ld); // S0star is now the choleksy decomposition of S0star^{-1};

		nustar = nclus_iter + nu0;


		ran_wish(nustar, S0star, nb, scr1, scr2, zerov, outrwish);

		cholesky(outrwish, nb, &ld);

		inverse_from_cholesky(outrwish, scr1, scr2, nb);
				
		for(b = 0; b < nb*nb; b++){
						
			Sig0_iter[b] = outrwish[b];
			Sig0I[b] = outrwish[b];
							
		}

						
		cholesky(Sig0I, nb, &ld);
					
		// I need to invert Astar to get an inverse wishart from the ran_wish function;
		inverse_from_cholesky(Sig0I, scr1, scr2, nb); 


		matrix_product(Sig0I, beta0_iter, Sig0Ibeta0, nb, 1, nb);



		//////////////////////////////////////////////////////////////////////////////////
		//
		// Update sig2 for each of the clusters
		//
		//////////////////////////////////////////////////////////////////////////////////
		
		for(k = 0; k < nclus_iter; k++){

//			Rprintf("k = %d =================== \n", k);
			sumsq = 0;
			for(j = 0; j < *nobs; j++){

				if(Si_iter[j] == k+1){
					xb = 0.0;
					for(b = 0; b < nb; b++){
						xb = xb + betah[k*(nb) + b]*X[j*(nb)+b];
					}
				
					sumsq = sumsq + (y[j] - xb)*(y[j] - xb);
				}

			}

//			Rprintf("sumsq = %f\n", sumsq);

			astar = 0.5*(nh[k]) + asig;
			bstar = 1/bsig + 0.5*sumsq;

//			Rprintf("astar = %f\n", astar);
//			Rprintf("bstar = %f\n", bstar);
						
			sig2h[k] = 1/rgamma(astar, 1/bstar);// 
		
		}

//		RprintVecAsMat("sig2h", sig2h, 1, nclus_iter);

//		RprintVecAsMat("sig2h", sig2h, 1, nclus_iter);

/*

		// Use a uniform prior
		for(k = 0; k < nclus_iter; k++){
//			Rprintf("k = %d =================== \n", k);
			so = sqrt(sig2h[k]);
			sn = rnorm(so,0.1);

//			Rprintf("so = %f\n", so);
//			Rprintf("sn = %f\n", sn);
//			Rprintf("sig2h[k] = %f\n", sig2h[k]);

			if(sn > A0 & sn < A){

				lln = 0.0, lln=0.0;
				for(j = 0; j < *nobs; j++){

					if(Si_iter[j] == k+1){
						xb = 0.0;
						for(b = 0; b < nb; b++){
							xb = xb + betah[k*(nb) + b]*X[j*(nb)+b];
						}
						
						llo = llo + dnorm(y[j], xb, so, 1);
						lln = lln + dnorm(y[j], xb, sn, 1);
					}
	
				}

//				Rprintf("llo = %f\n", llo);
//				Rprintf("lln = %f\n", lln);

				llo = llo + dunif(so, A0, A, 1);
				lln = lln + dunif(sn, A0, A, 1);

//				Rprintf("llo = %f\n", llo);
//				Rprintf("lln = %f\n", lln);

				llr = lln - llo;

//				Rprintf("llr = %f\n", llr);
				
				uu = runif(0,1);
				
				if(llr > log(uu)) sig2h[k] = sn*sn;
		
//				Rprintf("sig2h = %f\n", sig2h[k]);
			}
		}

*/
//		RprintVecAsMat("sig2h", sig2h, 1, nclus_iter);


		//////////////////////////////////////////////////////////////////////////////////
		//																				//
		// udpate M DP scale parameter;													//
		//																				//
		//////////////////////////////////////////////////////////////////////////////////
		if(*Mval == 0){
		
		
		}

		//////////////////////////////////////////////////////////////////////////////////
		//
		// evaluating likelihood that will be used to calculate LPML? 
		// (see page 81 Christensen Hansen and Johnson) 
		//
		//////////////////////////////////////////////////////////////////////////////////

		//
		//	Note the syntax for the pnorm function
		//	double pnorm(double x, double mu, double sigma, int lower_tail,int give_log);
		//

		if((i > (*burn-1)) & ((i+1) % (*thin) == 0)){

//			Rprintf("nclus = %d\n", nclus_iter);
//			RprintVecAsMat("muh", muh, 1, nclus_iter);
//			RprintVecAsMat("sig2h", sig2h, 1, nclus_iter);
//			RprintVecAsMat("y", y, 1, *nobs);


			for(j = 0; j < *nobs; j++){
//				Rprintf("j = %d\n", j);	

				s2tmp = sig2h[Si_iter[j]-1];

				xb = 0.0;
				for(b = 0; b < nb; b++){
					xb = xb + betah[(Si_iter[j]-1)*(nb) + b]*X[j*(nb)+b];
				}


				like_iter[j] = 	dnorm(y[j], xb, sqrt(s2tmp), 1);


//				Rprintf("like_iter = %f\n", like_iter[j]);
//				Rprintf("like_iter = %40.9f\n", like_iter[j]);

				// These are needed for WAIC
				mnlike[j] = mnlike[j] + exp(like_iter[j])/(double) nout; 
				mnllike[j] = mnllike[j] + (like_iter[j])/(double) nout;

					
				CPOinv[j] = CPOinv[j] + (1/(double) nout)*
					                  		(1/exp(like_iter[j]));

//				Rprintf("CPO = %f\n", CPO[j]);

			}

						
		}



	
		//////////////////////////////////////////////////////////////////////////////////
		//																				//
		// Save MCMC iterates															//
		//																				//
		//////////////////////////////////////////////////////////////////////////////////
		if((i > (*burn-1)) & ((i+1) % *thin ==0)){

//			Rprintf("ii = %d\n", ii);

			nclus[ii] = nclus_iter;
//			Rprintf("nclus = %d\n", nclus[ii]); 


			for(b = 0; b < nb; b++){
				beta0[ii*(nb) + b] = beta0_iter[b];
			}

//			RprintVecAsMat("beta0", beta0, nout, nb);
				
			for(b = 0; b < nb*nb; b++){

				Sig0[ii*(nb*nb) + b] = Sig0_iter[b];

			}



			if(*Mval == 0){
				M[ii] = M_iter;
			}
		

			for(j = 0; j < *nobs; j ++){
				xb = 0.0;
				for(b = 0; b < nb; b++){
					beta[(ii*(*nobs) + j)*nb + b] = betah[(Si_iter[j]-1)*nb + b];
										
					xb = xb + betah[(Si_iter[j]-1)*nb + b]*X[j*(nb)+b];


				}

				sig2[ii*(*nobs) + j] = sig2h[Si_iter[j]-1];
				Si[ii*(*nobs) + j] = Si_iter[j];
				llike[ii*(*nobs) + j] = like_iter[j];
				fitted[ii*(*nobs) + j] = xb;
				
			}

/*			// Store density values
			for(b = 0; b < nb*nb; b++) scr2[b] = Sig0_iter[b];
			cholesky(scr2, nb , &ld);
			ran_mvnorm(beta0_iter, scr2, nb, scr1, betadraw);
			sig2draw = 1/rgamma(asig, bsig);

			for(jj = 0; jj < *nzgrid; jj++){
				for(g = 0; g < *nxpred; g++){
					dval = 0.0;				
					for(k = 0; k < nclus_iter; k++){
						xb = 0.0;
						for(b = 0; b < nb; b++){
							xb = xb + betah[k*nb + b]*Xpred[g*nb+b];
						}
					
						dval = dval + ((double) nh[k])/(*nobs + M_iter)*
						              dnorm(zgrid[jj], xb, sqrt(sig2h[k]),0);

					  	cumpval = cumpval + ((double) nh[k])/(*nobs + M_iter)*
					    			pnorm(zgrid[jj], xb, sqrt(sig2h[k]),1,0);
					}
			

					xb = 0.0;
					for(b = 0; b < nb; b++){
						xb = xb + betadraw[b]*Xpred[g*nb + b];
					}

					dval = dval + (M_iter)/(*nobs + M_iter)*
					             dnorm(zgrid[jj], xb, sqrt(sig2draw),0);

					cumpval = cumpval + (M_iter)/(*nobs + M_iter)*
				  				pnorm(zgrid[jj], xb, sqrt(sig2draw),1,0);

					density[(ii*(*nzgrid) + jj)*(*nxpred) + g] = dval;
					cumdensity[(ii*(*nzgrid) + jj)*(*nxpred) + g] = cumpval;

				
				}

			}
*/
			ii = ii+1;
		
		}



	}



	lpml_iter=0.0;
	for(j = 0; j < *nobs; j++){
//		Rprintf("j = %d\n", j);

//		Rprintf("CPO = %f\n", CPO[j*(*ntime)+t]);
				
		lpml_iter = lpml_iter - log(CPOinv[j]);
				
	}
//	Rprintf("nout_0 = %d\n", nout_0);
	lpml[0] = lpml_iter;
//	Rprintf("lpml_iter = %f\n", lpml_iter);


	elppdWAIC = 0.0;
//	RprintVecAsMat("mnlike",mnlike, 1,(*nsubject)*(ntime1));
//	RprintVecAsMat("mnllike", mnllike, 1, (*nsubject)*(ntime1));
		
	for(j = 0; j < *nobs; j++){
		elppdWAIC = elppdWAIC + (2*mnllike[j] - log(mnlike[j]));  
	}
	waic[0] = -2*elppdWAIC; 

	PutRNGstate();
	



}


