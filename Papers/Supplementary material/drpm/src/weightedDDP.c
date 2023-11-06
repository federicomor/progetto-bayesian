/*************************************************************
 * Copyright (c) 2009 Garritt Leland Page
 *
 * This file contains C code for an MCMC algorithm 
 * associataed with the the DP prior lab mean model.  This will
 * be used to perfrom the simulation study for the paper
 * that I worked on with David Dunson.
 * An outline of the model follows
 *
 *        y_{ij} ~ N_p(mu_i, Sigma_i)
 *          mu_i ~ P
 *             P ~ DP(alpha, P0)
 *            P0 ~ N_p(mu0, Sigma0)
 *		     mu0 ~ N_p(m, S)
 *                    
 *
 * Inputs of the function are
 * draws - number of mcmc iterates
 * burn - number of initial draws to remove
 * thin - indicates which draws to be retained
 *
 * nobs - number of observations (rows in Ymat)
 * ndim - number of variables (columns in Ymat)
 * N - upper bound on the number of components in the infinite mixture
 * Ymat - nobs x ndim matrix of observations.  the first columns will 
          correspond to the response and the subsequent columns correspond 
          to covariates.
 * 
 * The following are scratch vectors to hold MCMC iterates
 * mu - (draws-burn)/thin x ndim x N array containing MCMC iterates for mu
 * Sigma - (draws-burn)/thin x 
 * 
 *************************************************************/
 
#include "matrix.h"
#include "Rutil.h"

#include <R_ext/Lapack.h>
#include <R.h>
#include <Rmath.h>


//#define size 3

void WDDP(int *draws, int *burn, int *thin, int *nobs, int *dim, int *N,
          double *Ymat, double *priorparms,
		  double *mu, double *Sigma, int* Si, int *nclus, 
		  double *mu0, double *Sigma0, double* sbw,
		  double *llike, double *fitted, double *lpml, double *waic){
	


    int i, j, t, tt, d, ii, T, h;
	int ncp=(*dim)*(*dim);
	int df;
	
//	RprintVecAsMat("Ymat", Ymat, *nobs, *dim);
//	Rprintf("N = %d\n", *N);
	int nout;
//	double matrix[dim*dim]  ;

    nout = (*draws - *burn)/(*thin);

	double uu, ld;

	T = (*nobs)*(*dim);//Total number of observations
	

	double alpha;
	
	alpha = 1.0;

//	RprintVecAsMat("SmuI = ", SmuI,dim,dim);
//	RprintVecAsMat("csig1 = ", csig1,1,dim);
//	RprintVecAsMat("csig2 = ", csig2,1,dim);

	// ===================================================================================
	// 
	// Memory vectors to hold MCMC iterates for cluster specific parameters
	// 
	// ===================================================================================
	int Si_iter[*nobs]; for(j=0; j<*nobs; j++) Si_iter[j]=rbinom(2 , 0.5) + 1;
//    RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);
	int nclus_iter; 

	double *mu_iter = R_VectorInit((*N)*(*dim), 0.0);
	double *mu0_iter = R_VectorInit((*N)*(*dim), 0.0);
	double *Sigma_iter = R_VectorInit((*N)*(ncp), 0.0);
	double *Sigma0_iter = R_VectorInit((*N)*(ncp), 0.0);


	// scratch vectors to use in the matrix functions found in matrix.c
	double *tmp1 = R_VectorInit(*dim, 0.0);
	double *tmp2 = R_VectorInit(*dim, 0.0);
	double *tmp3 = R_VectorInit(*dim, 0.0);
	double *zers = R_VectorInit(*dim, 0.0);

	// Stuff that I need for updating Sigma_j;
	double *sumy_muy_mup = R_Vector(ncp);
	double *Astar = R_Vector(ncp);	
	double *outrwish = R_Vector(ncp);
	double *iSig = R_VectorInit((*N)*(ncp), 0.0);

	// Stuff for updating mu;
	double *sumy = R_Vector(*dim);
	double *Mstar = R_Vector(*dim);
	double *Sstar = R_Vector(*dim);
	double *outrmvnorm = R_Vector(*dim);
	double *iSig0 = R_VectorInit(ncp, 0.0);
    for(t=0;t<*dim;t++){for(tt=0;tt<*dim;tt++){if(t==tt){iSig0[t*(*dim)+tt]=1.0;}}}

	// Stuff for updating Si
	int nh;
	double maxp, den, cprobh; 
	double ldet[*N], ph[*N];
	double *sbweight = R_VectorInit(*N, 1/(*N));
	double *mutmp = R_Vector(*dim);
	double *ytmp = R_Vector(*dim);
	double *iSigtmp = R_Vector(ncp);


	// Stuff for updating the stick-breaking weights
	double ccount, cumccount, csbw;
	double astar, bstar;

	// Stuff that  I need for updating mu0;
	double *summu = R_Vector(*dim);

	// Stuff that I need for updating V;
	double* summu_mu0 = R_VectorInit(ncp, 0.0);
	
		


	// Stuff to compute lpml, likelihood, WAIC, and Rao-Blackwellized density values
	// and other model evaluation metrics
	double lpml_iter, elppdWAIC, muy, sxdens;
	double Sigxy[(*dim)*(*dim)],  Sigxx[(*dim)*(*dim)], mux[*dim], xvec[*dim];
	double condmn[*N], condvar[*N], xdens[*N], condweight[*N];
	double *CPOinv = R_VectorInit(*nobs, 0.0);
	double *like_iter = R_VectorInit(*nobs, 0.0);
	double *mnlike = R_VectorInit(*nobs, 0.0);
	double *mnllike = R_VectorInit(*nobs, 0.0);
    double *ispred = R_VectorInit(*nobs, 0.0);


	// Stuff for priors

	// priors for mu0
	double *M0=R_VectorInit(*dim, priorparms[0]);
	double *S0=R_VectorInit(ncp, 0);
	double *iS0=R_VectorInit(ncp, 0);

	for(t = 0; t < *dim; t++){
	  for(tt = 0; tt < *dim; tt++){
	    if(t == tt){
	      S0[t*(*dim) + tt] = priorparms[1];
	      iS0[t*(*dim) + tt] = (1/priorparms[1]);
        }
	  }
    }
//    RprintVecAsMat("M0", M0, 1, *dim);
//    RprintVecAsMat("S0", S0, *dim, *dim);
//    RprintVecAsMat("iS0", iS0, *dim, *dim);
    
	// Need the prior matrix A for distributions of Sigma and Sigma0;
	double* A = R_VectorInit(ncp,0.0);
	for(d = 0; d < ncp; d++) if(d % (*dim+1) == 0.0) A[d] = priorparms[3]; 
//	RprintVecAsMat("A =", A, *dim, *dim);
	
	double nu1 = priorparms[2];
//	Rprintf("nu1 = %f\n", nu1);	


    Rprintf("Prior values: M0 = %.2f, S0 = %.2f\n \t A = %.2f, nu = %.2f\n",
      	           priorparms[0], priorparms[1], priorparms[3], nu1);
	
	GetRNGstate();

    ii = 0;

	// start of the mcmc algorithm;
    for(i = 0; i < *draws; i++){

      if(i % 1 == 10000) Rprintf("mcmc iter = %d\n", i);

      /************************************************************************** 
       *                                                                        *
       *  update Sigma_j (component specific covariance matrix)                 *
       *                                                                        *
       **************************************************************************/ 

	  for(h = 0; h < *N; h++){
//		Rprintf("h = %d\n", h);

        for(d = 0; d < ncp;  d++) sumy_muy_mup[d] = 0.0;
        
        nh = 0;			
        for(j = 0; j < *nobs; j++){
        
          if(Si_iter[j] == h+1){
            nh = nh+1;
            for(t = 0; t < *dim; t++){
              for(tt = 0; tt < *dim; tt++){
                sumy_muy_mup[t*(*dim) + tt] = sumy_muy_mup[t*(*dim) + tt] +  
                     (Ymat[j*(*dim) + t] - mu_iter[h*(*dim) + t])*
                     (Ymat[j*(*dim) + tt] - mu_iter[h*(*dim) + tt]);
        	  }			
        	}
          }
        }
        
//        RprintVecAsMat("sumy_muy_mup", sumy_muy_mup, *dim, *dim);
//        Rprintf("nh = %d\n", nh);
        for(d = 0; d < ncp; d++){
          Astar[d] = sumy_muy_mup[d] + A[d];
        }
        			
        cholesky(Astar, *dim, &ld);
        
//        Rprintf("ld = %f\n", ld);
        // I need to invert Astar to get an inverse wishart from the ran_wish function;
        inverse_from_cholesky(Astar, tmp1, tmp2, *dim); 
        			
        cholesky(Astar, *dim, &ld); // Astar is now the choleksy decomposition of Astar^{-1};
        
        df = nh + nu1;
        
        ran_wish(df, Astar, *dim, tmp1, tmp2, zers, outrwish);
        			
        for(d = 0; d < ncp; d++){
          iSig[h*ncp + d] = outrwish[d];
        }	
//        RprintVecAsMat("outrwish", outrwish, *dim, *dim);
        cholesky(outrwish, *dim, &ld);
        ldet[h] = -ld; 
        // Note that I use -ld rather than ld since I need log det(Sigma) and the 
        // cholesky is being applied to Sigma^{-1} and hence produces log det(Sigma^{-1})

//        Rprintf("ldet[h] = %f\n",ldet[h]);

        // When I input the matrix outrwish into R and use the solve function
        // the results are not exactly the same as those I get from using the inverse_from_choleksy 
        // function here.  This might be something to be considered about.  Although the cholesky is
        // pretty much right on.  Also Astar had a few minor issues but I attribute it to round off error;
        
        inverse_from_cholesky(outrwish, tmp1, tmp2, *dim);
//        RprintVecAsMat("outrwish", outrwish, *dim, *dim);
        	
        for(d = 0; d < ncp; d++){
          Sigma_iter[h*ncp + d] = outrwish[d];
        }
      }
//      RprintVecAsMat("Sigma", Sigma_iter, *N, ncp);

      /************************************************************************** 
       *                                                                        *
       *  update mu_j (component specific mean vector)                          *
       *                                                                        *
       **************************************************************************/ 
      for(t = 0; t < *dim; t++) summu[t]=0.0; 

      for(h = 0; h < *N; h++){

//        Rprintf("h = %d\n", h);
        
        for(t = 0; t < *dim; t++) sumy[t] = 0.0;
//        Rprintf("h+1 = %d\n", h+1);
//        RprintVecAsMat("sumytmp", sumytmp, 1, dim);
        nh = 0;
        for(j = 0; j < *nobs; j++){
          if(Si_iter[j] == h+1){
        	nh = nh + 1;		
            for(t = 0; t < *dim; t++){
        	  sumy[t] = sumy[t] + Ymat[j*(*dim)+t];
        	}
          }
        }
//        RprintVecAsMat("sumy", sumy, 1, *dim);

		for(d = 0; d < ncp; d++){
		  iSigtmp[d] = iSig[h*ncp + d];
		  Sstar[d] = ((double) nh)*iSigtmp[d] + iSig0[d]; // recall that Sigma and Tau are inverses.;
		}


        
		cholesky(Sstar, *dim, &ld);
		inverse_from_cholesky(Sstar, tmp1, tmp2, *dim);
//		RprintVecAsMat("Sstar", Sstar, dim, dim);
        
        
        matrix_product(iSig0, mu0_iter, tmp1, *dim, 1, *dim);
        
        matrix_product(iSigtmp, sumy, tmp2, *dim, 1, *dim); 				

		for(t = 0; t < *dim; t++){
		  tmp3[t] = tmp1[t] + tmp2[t];
		}

		matrix_product(Sstar, tmp3, Mstar, *dim, 1, *dim); //recall Sstar is an inverse

		// to use the ran_mvnorm function I need to get the cholesky decomposition of the 
		// covariance matrix.  In this case it is Sstar (which recall has already been inverted so 
		// it is the correct covariance matrix) 
		cholesky(Sstar, *dim , &ld);
																
		ran_mvnorm(Mstar, Sstar, *dim, tmp1, outrmvnorm);
					
		for(t = 0; t < *dim; t++){
		  mu_iter[h*(*dim) + t] = outrmvnorm[t];
		  summu[t] = summu[t] + mu_iter[h*(*dim) + t];
		}
		
	  }
//	  RprintVecAsMat("mu", mu_iter, *N, *dim);										


      /************************************************************************** 
       *                                                                        *
       *  Next I update mu0.  This is the mean of the baseline distribution     *
       *  found in the DP prior. The complete conditional is a MVN normal       *    
       *  distribution.                                                         *
       *                                                                        *
       **************************************************************************/ 

	  for(d = 0; d < ncp; d++){
	    Sstar[d] = (*N)*iSig0[d] + iS0[d]; 
	  }
				
	  cholesky(Sstar, *dim, &ld);
//    RprintVecAsMat("SstarChol = ", Sstar,dim,dim);

	  inverse_from_cholesky(Sstar, tmp1, tmp2, *dim); // Sstar is now an inverse 
//	  RprintVecAsMat("Sstar", Sstar,*dim,*dim);

 //   RprintVecAsMat("summu = ", summu,1,*dim);

	  matrix_product(iSig0, summu, tmp1, *dim, 1, *dim);
	  matrix_product(iS0, M0, tmp2, *dim, 1, *dim);
//    RprintVecAsMat("SiM = ", SiM,1,dim);
			
	  for(t = 0; t < *dim; t++){
	    tmp3[t] = tmp1[t] + tmp2[t];
	  }
				
	  matrix_product(Sstar, tmp3, Mstar, *dim, 1, *dim);  // recall that Sstar is an inverse
//    RprintVecAsMat("Mstar = ", Mstar,1,dim);


      // This cholesky of Sstar doesn't match what I get when using chol in R
	  cholesky(Sstar, *dim, &ld); // I need the cholesky decomposition of the cov mat to use ran_mvnorm; 
//    RprintVecAsMat("Sstar = ", Sstar, dim ,dim);

	  ran_mvnorm(Mstar, Sstar, *dim, tmp1, outrmvnorm);

	  for(t = 0; t < *dim; t++){
	    mu0_iter[t] = outrmvnorm[t];
	  }
//	  RprintVecAsMat("mu0", outrmvnorm, 1, *dim);										


      /************************************************************************** 
       *                                                                        *
       *  Next I update Sig0.   The cc is IW.                                   *
       *																		*      
       **************************************************************************/ 

      for(d = 0; d < ncp; d++) summu_mu0[d] = 0;

      for(h = 0; h < *N; h++){
        for(t = 0; t < *dim; t++){
          for(tt = 0; tt < *dim; tt++){
            summu_mu0[t*(*dim) + tt] = summu_mu0[t*(*dim) + tt] +  
                     (mu0_iter[t] - M0[t])*
                     (mu0_iter[tt] - M0[tt]);
          }			
        }
      }

	  for(d = 0; d < ncp; d++){
	    Astar[d] = summu_mu0[d] + A[d];
	  }
			
	  cholesky(Astar, *dim, &ld);

	  inverse_from_cholesky(Astar, tmp1, tmp2, *dim);

	  cholesky(Astar, *dim, &ld);


	  df = (*N) + (*dim + 1);

	  ran_wish(df, Astar, *dim, tmp1, tmp2, zers, outrwish);

      for(d = 0; d < ncp; d++){
        iSig0[d] = outrwish[d];
      }

	  cholesky(outrwish, *dim, &ld);
	  inverse_from_cholesky(outrwish, tmp1, tmp2, *dim);

	  for(d = 0; d < ncp; d++){
						
	    Sigma0_iter[d] = outrwish[d];
							
	  }

//      RprintVecAsMat("Sigma0 = ", Sigma0_iter,*dim,*dim); 
			


       
      /************************************************************************** 
       *                                                                        *
       *  Next I updated the allocation of units to components.                 * 
       *  Recall that I am using the block Gibbs sampler which                  *
       *  requires setting an upper bound to the number of atoms.               *
       *  The allocation parameter is Si                                        *
       *                                                                        *
       **************************************************************************/ 
      
	  // Now I need to evaluate the MVN normal density for each of the component
	  // specific atoms
			
      for(j = 0; j < *nobs; j++){

//        Rprintf("j = %d\n", j);
        for(h = 0; h < *N; h++){
//          Rprintf("h = %d\n", h);
          for(t = 0; t < *dim; t++){
										
		    ytmp[t] = Ymat[j*(*dim) + t];
			mutmp[t] = mu_iter[h*(*dim) + t];
		  }
		  
		  for(d = 0; d < ncp; d++){
		    iSigtmp[d] = iSig[h*ncp + d];
		  }
//          RprintVecAsMat("ytmp", ytmp, 1, *dim);
//          RprintVecAsMat("mutmp", mutmp, 1, *dim);
//          RprintVecAsMat("iSigtmp", iSigtmp, *dim, *dim);
//          Rprintf("dmvnorm(ytmp, mutmp, iSigtmp, *dim, ldet[h], tmp1, 1) = %f\n", dmvnorm(ytmp, mutmp, iSigtmp, *dim, ldet[h], tmp1, 1));

		  ph[h] = dmvnorm(ytmp, mutmp, iSigtmp, *dim, ldet[h], tmp1, 1) + log(sbweight[h]);
		  
		}								

//        RprintVecAsMat("ph", ph, 1, *N);
        maxp = ph[0];
        for(h = 1; h < *N; h++) if(ph[h] > maxp) maxp = ph[h];

        den = 0.0;
		for(h = 0; h < *N; h++){
		  ph[h] = exp(ph[h] - maxp);
		  den = den + ph[h];
		}


		uu = runif(0.0,1.0);
		cprobh= 0.0;;
		for(h = 0; h < *N; h++){
		  cprobh = cprobh + ph[h]/den;
		  if(uu < cprobh){
		    Si_iter[j] = h+1;
		    break;	
		  }
		}
      }
      
      nclus_iter = 0;
      for(j = 0; j < *nobs; j++){
        for(h = 0; h < *N; h++){
          if(Si_iter[j] == h+1){
            nclus_iter = nclus_iter + 1;
            break;
          }
        }
      }
//      RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);
      
      /************************************************************************** 
       *                                                                        *
       *  Next I update the stick-breaking weights.  These are draws from a     *
       *  beta.  We draw N-1 of them and then make the Nth weight 1.            *
       *                                                                        *
       **************************************************************************/ 
      csbw = 1.0;
	  for(h = 0; h < *N; h++){
	    ccount = 0.0;
	    cumccount = 0.0;
		for(j = 0; j < *nobs; j++){

//		  Rprintf("Si_iter[j] = %d\n", Si_iter[j]);
						
		  if(Si_iter[j] == h+1){
		    ccount = ccount + 1.0;
		  }

		  if(Si_iter[j] > h+1){
		    cumccount = cumccount + 1.0;
		  }
		}

//		Rprintf("ccount = %f\n", ccount);
//		Rprintf("cumccount = %f\n", cumccount);
		astar = 1.0 + ccount;
		bstar = alpha + cumccount;
					
//		Rprintf("astar = %f\n", astar);
//		Rprintf("bstar = %f\n", bstar);
		sbw[h] = rbeta(astar,bstar);

		if(h == *N-1) sbw[h] = 1.0;

		if(h == 0){
		  sbweight[h] = sbw[h];
		}
		if(h > 0){
		  csbw = csbw*(1-sbw[h-1]);
		  sbweight[h] = sbw[h]*csbw;
		}			
	  }

//       RprintVecAsMat("sbw", sbw, 1, *N);
//       RprintVecAsMat("sbweight", sbweight, 1, *N);

    ////////////////////////////////////////////////////////////////////////////////////////////
    //
    // in sample prediction to assess model fit
    // Use the conditional specification of the model
    //
    ////////////////////////////////////////////////////////////////////////////////////////////
    if((i > (*burn-1)) & (i % (*thin) == 0)){
      
      for(j = 0; j < *nobs; j++){
//        Rprintf("j = %d\n", j);
        // find the sum of marginal density values
        for(t = 1; t < *dim; t++){
          xvec[t-1] = Ymat[j*(*dim) + t];
        }
//        RprintVecAsMat("xvec", xvec, 1, *dim-1);
        
        sxdens = 0.0;
        for(h = 0; h < *N; h++){
//          Rprintf("h = %d\n", h);
          muy = mu_iter[h*(*dim) + 0];
          for(t = 1; t < *dim; t++){
            mux[t-1] = mu_iter[h*(*dim) + t];
            for(tt = 1; tt < *dim ; tt++){
//              Rprintf("h*ncp + t*(*dim) + tt = %d\n", h*ncp + t*(*dim) + tt);
              Sigxx[(t-1)*(*dim - 1) + (tt-1)] = Sigma_iter[h*ncp + t*(*dim) + tt];
            }
            for(tt = 0; tt < *dim-1; tt++){
              Sigxy[tt] = Sigma_iter[h*ncp + t*(*dim) + tt];
            }
          }

//          RprintVecAsMat("Sigxx", Sigxx, *dim-1, *dim-1);
//          RprintVecAsMat("Sigxy", Sigxy, 1, *dim-1);

	      cholesky(Sigxx, *dim-1, &ld);
	      inverse_from_cholesky(Sigxx, tmp1, tmp2, *dim-1);
          
          xdens[h] = dmvnorm(xvec, mux, Sigxx, *dim-1, ld, tmp1, 0)*(sbweight[h]);
          sxdens = sxdens + xdens[h];

          condmn[h] = 0.0;
          condvar[h] = 0.0;
          for(t = 0; t < *dim - 1; t++){
            for(tt = 0; tt < *dim - 1; tt++){
              condmn[h] = condmn[h] + Sigxy[t]*Sigxx[t*(*dim-1)+tt]*(xvec[tt] - mux[tt]); 
              condvar[h] = condvar[h] + Sigxy[t]*Sigxx[t*(*dim-1)+tt]*Sigxy[tt];
            }
          }
//          Rprintf("condmn = %f\n", condmn[h]);
//          Rprintf("condvar = %f\n", condvar[h]);
          condmn[h] = condmn[h] + muy;
          condvar[h] =  Sigma_iter[h*ncp + 0] - condvar[h];

        }
        
//        RprintVecAsMat("xdens", xdens, 1, *N);
//        RprintVecAsMat("condmn", condmn, 1, *N);
//        RprintVecAsMat("condvar", condvar, 1, *N);

        ispred[j] = 0.0;
        for(h = 0; h < *N; h++){
          condweight[h] = xdens[h]/sxdens;
          ispred[j] = ispred[j] + condweight[h]*condmn[h];
        }
//        ispred[j] = condmn[Si_iter[j] - 1];
      
      
        /////////////////////////////////////////////
        //
        // Compute the CPO and lpml using the mixture
        //
        /////////////////////////////////////////////
        like_iter[j] = 0.0;
        for(h = 0; h < *N; h++){
		  like_iter[j] = like_iter[j] + 
		                  condweight[h]*
		                  dnorm(Ymat[j*(*dim) + 0], condmn[h], sqrt(condvar[h]), 0);
        }

        // These are needed for WAIC
        mnlike[j] = mnlike[j] + (like_iter[j])/(double) nout;
        mnllike[j] = mnllike[j] + log(like_iter[j])/(double) nout;

        CPOinv[j] = CPOinv[j] + (1/(double) nout)*(1/like_iter[j]);

     }
      
    }


/*
	// ===================================================================================
	//
	// Out of sample prediction
	//
	// ===================================================================================

	if(i > (*burn-1) & i % (*thin) == 0){



//	  RprintVecAsMat("sbweight", sbweight, 1, *N); 
//	  RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs); 
	  llike = 0.0;	
	  for(pp = 0; pp < *npred; pp++){

//	    Rprintf("pp =============== %d\n", pp);
							
	    for(t = 0; t < *dim; t++) xpredobstmp[t] = xpredobs[pp*(*dim) + t];

//	    RprintVecAsMat("xpredobstmp", xpredobstmp, 1, *m); 

	    cwdensMix=0.0;
	    for(h = 0; h < *N; h++){
//	      Rprintf("h = %d\n", h);
	      for(t = 0; t < *dim; t++){
		    ytmp[t] = Ymat[j*(*dim) + t];
		    mutmp[t] = mu_iter[h*(*dim) + t];
		  }
		  for(d = 0; d < ncp; d++){
		    iSigtmp[d] = iSigma[h*ncp + d]
		  }

//		  RprintVecAsMat("mutmp = ", mutmp, 1, *dim);
//		  RprintVecAsMat("iSig = ", iSigtmp, *dim, *dim);
//		  Rprintf("det = %f\n", det);
//		  Rprintf("dmvnorm = %f\n",dmvnorm(ytmp, mutmp, iSigtmp, *dim, ldet[h], tmp1, 1));
//		  Rprintf("sbweight = %f\n", sbweight[h]);								

		  mvnDens[h] = dmvnorm(ytmp, mutmp, iSigtmp, *dim, ldet[h], tmp1, 1) + log(sbweight[h]);

		  llike = llike + mvnDens[h];
	    }


		for(h = 0; h < *N; h++) snumSi[h] = mvnDens[h];
					
				R_rsort (snumSi,  *N) ;
					
//				RprintVecAsMat("snumSi ", snumSi, 1, *N);
					
				maxnumSi = snumSi[*N-1];
					
//				Rprintf("maxnumSi = %f\n", maxnumSi);
					
//				RprintVecAsMat("numSi", numSi, 1, *N);

				for(h = 0; h < *N; h++){
						
					mvnDens[h] = exp(mvnDens[h] - maxnumSi);
//					numSi[h] = pow(exp(numSi[h] - maxnumSi), (1 - exp(-1*(i+1))));
					cwdensMix = cwdensMix + mvnDens[h];
					
				}

//				RprintVecAsMat("mvnDens", mvnDens, 1, *N);
//				Rprintf("cwdensMix = %f\n", cwdensMix);
			
				for(h = 0; h < *N; h++){
				
				  condwght[h] = mvnDens[h]/cwdensMix;
					
				}

//				RprintVecAsMat("sbweight", sbweight, 1, *N); 
//				RprintVecAsMat("condwght", condwght, 1, *N);
//				Rprintf("ypredobs = %f\n", ypredobs[pp]);
//				Rprintf("tau2ITER = %f\n", tau2ITER);
			
								
//				RprintVecAsMat("condwght", condwght, 1, *N);
//				RprintVecAsMat("thetayITER", thetayITER, *N, *C);

				// Now for predictions which we can treat in two ways.  one is the
				// expected value of the other is from the "posterior predictive". First
				// the expected value and then the posterior predictive

				for(c = 0; c < *C; c++) ypprob[c] = 0.0;;

				for(c = 0; c < *C; c++){
//					Rprintf("h = %d\n", h);

					for(h = 0; h < *N; h++){

						ypprob[c] = ypprob[c] + condwght[h]*thetayITER[h*(*C) + c];
					}
					
				}		

		
				for(c = 0; c < *C; c++) fprintf(ypredprobMCMC, "%10.9f\n", ypprob[c]);


//				RprintVecAsMat("ypprob", ypprob, 1, *C);
							
				uu = runif(0.0,1.0);

				cprobSi = 0.0;;
				for(c = 0; c < *C; c++){
//					Rprintf("h = %d\n", h);

					cprobSi = cprobSi + ypprob[c];

//					Rprintf("cprobSi = %f\n", cprobSi);

					if (uu < cprobSi){
							
						ypredITER[pp] = c+1;
						break;	
					}
					
				}		
							
//				Rprintf("cprobSi = %f\n", cprobSi);
//				Rprintf("ypred = %f\n", ypredITER[pp]);											

				fprintf(ypredMCMC,"%10.9f\n", ypredITER[pp]); 
								
			}


			// =====================================================================================
			//
			// predict in sample ys.  
			// compute the cluster weights
			// compute the log pseudo marginal likelihood
			//
			// =====================================================================================
			
			
//			RprintVecAsMat("ypred", ypredITER, 1, *npred);				


			// obtain the conditional weights to fit the model to the cloud of points

			for(kk = 0; kk < (ssdimITER)*(ssdimITER); kk++){
				sig2mat[kk] = 0.0;

				if(kk % (ssdimITER + 1) == 0) sig2mat[kk] = 1/sig2ITER;

			}

//			RprintVecAsMat("sig2mat", sig2mat, ssdimITER, ssdimITER);

			det = pow(sig2ITER, ssdimITER);

//			Rprintf("det = %f\n", det);

//			ff = factorial(10);
//			Rprintf("ff = %d\n", ff);

			for(kk = 0; kk < (*m)*(*m); kk++){
				sig2mat2[kk] = 0.0;

				if(kk % (*m + 1) == 0) sig2mat2[kk] = 1/sig2ITER;

			}

//			RprintVecAsMat("sig2mat", sig2mat, ssdimITER, ssdimITER);

			det2 = pow(sig2ITER, *m);


			// =====================================================================================
			//
			// predict in sample ys.  
			//
			// =====================================================================================

//			llike = 0.0;			
			for(j = 0; j < *n; j++){

//				Rprintf("j = %d\n",j);

				for(k = 0; k < *m; k ++) xtmp[k] = xmat[j][k];

//				RprintVecAsMat("xtmp", xtmp, 1, *m);

				matrix_product(tUmatITER, xtmp, tUkx, ssdimITER, 1, *m);

//				RprintVecAsMat("tUkxpred", tUkxpred, 1, ssdimITER); 


				cwdensMix=0.0;
				for(h = 0; h < *N; h++){
						
//					Rprintf("h = %d\n", h);
					for(k = 0; k < ssdimITER; k++){
								
						mutmp[k] = thetaxITER[h*(ssdimITER) + k];
									
					}

//					RprintVecAsMat("tUkxpred = ", tUkxpred, 1, ssdimITER);
//					RprintVecAsMat("mutmp = ", mutmp, 1, ssdimITER);
//					RprintVecAsMat("sig2mat = ", sig2mat, ssdimITER, ssdimITER);
//					Rprintf("det = %f\n", det);
//					Rprintf("dmvnorm = %f\n",dmvnorm(tUkxpred, mutmp, sig2mat, ssdimITER, log(det), scratch1, 0));
//					Rprintf("sbweight = %f\n", sbweight[h]);								

					mvnDens[h] = log(sbweight[h]) + dmvnorm(tUkx, mutmp, sig2mat, ssdimITER, log(det), scratch1, 1);
//					mvnDens[h] = sbweight[h]*dmvnorm(tUkxpred, mutmp, sig2mat, ssdimITER, log(det), scratch1, 0);

					//cwdensMix = cwdensMix + mvnDens[h];
					llike = llike + mvnDens[h];
				}


				for(h = 0; h < *N; h++) snumSi[h] = mvnDens[h];
					
				R_rsort (snumSi,  *N) ;
					
//				RprintVecAsMat("snumSi ", snumSi, 1, *N);
					
				maxnumSi = snumSi[*N-1];
					
//				Rprintf("maxnumSi = %f\n", maxnumSi);
					
//				RprintVecAsMat("numSi", numSi, 1, *N);

				for(h = 0; h < *N; h++){
						
					mvnDens[h] = exp(mvnDens[h] - maxnumSi);
//					numSi[h] = pow(exp(numSi[h] - maxnumSi), (1 - exp(-1*(i+1))));
					cwdensMix = cwdensMix + mvnDens[h];
					
				}

//				RprintVecAsMat("mvnDens", mvnDens, 1, *N);
//				Rprintf("cwdensMix = %f\n", cwdensMix);
			
				for(h = 0; h < *N; h++){
				
					condwght[h] = mvnDens[h]/cwdensMix;
					
				}

//				RprintVecAsMat("sbweight", sbweight, 1, *N); 
//				RprintVecAsMat("condwght", condwght, 1, *N);
//				Rprintf("ypredobs = %f\n", ypredobs[pp]);
//				Rprintf("tau2ITER = %f\n", tau2ITER);
			
								
//				RprintVecAsMat("condwght", condwght, 1, *N);
//				RprintVecAsMat("thetayITER", thetayITER, *N, *C);

				// Now for predictions which we can treat in two ways.  one is the
				// expected value of the other is from the "posterior predictive". First
				// the expected value and then the posterior predictive

				for(c = 0; c < *C; c++) ypprob[c] = 0.0;;

				for(c = 0; c < *C; c++){
//					Rprintf("h = %d\n", h);

					for(h = 0; h < *N; h++){

						ypprob[c] = ypprob[c] + condwght[h]*thetayITER[h*(*C) + c];
					}
					
				}		


//				RprintVecAsMat("ypprob", ypprob, 1, *C);
							
				uu = runif(0.0,1.0);

				cprobSi = 0.0;;
				for(c = 0; c < *C; c++){
//					Rprintf("h = %d\n", h);

					cprobSi = cprobSi + ypprob[c];

//					Rprintf("cprobSi = %f\n", cprobSi);

					if (uu < cprobSi){
							
						yinsamplepredITER[j] = c+1;
						break;	
					}
					
				}		
							
//				Rprintf("cprobSi = %f\n", cprobSi);
//				Rprintf("ypred = %f\n", ypredITER[pp]);											

				fprintf(yinsamplepredMCMC,"%10.9f\n", yinsamplepredITER[j]); 


				
			// =====================================================================================
			//
			// log pseudo likelihood and individual weights in conditional mixture
			//
			// =====================================================================================


				for(c = 0; c < *C; c++){ 
					ytmp[c] = 0;
					if(ys[j] == c+1) ytmp[c] = 1;
				}
				for(k = 0; k < *m; k++) xtmp[k] = xmat[j][k];

//				RprintIVecAsMat("ytmp", ytmp, 1, *C);

				matrix_product(tUmatITER, xtmp, tUkx, ssdimITER, 1, *m);

//				RprintVecAsMat("tUkx", tUkx, 1, ssdimITER);

						
				cwdensMix=0.0;
				loglike = 0.0;
//				Rprintf("*N = %d\n", *N);

				for(h = 0; h < *N; h++){
						
//					Rprintf("h = %d\n", h);

					for(k = 0; k < ssdimITER; k++){
								
						mutmp[k] = thetaxITER[h*(ssdimITER) + k];
					}
					
					for(c = 0; c < *C; c++){

						nutmp[c] = thetayITER[h*(*C) + c];

					}

					matrix_product(UmatITER, mutmp, Umu, *m, 1, ssdimITER);
				
//					RprintVecAsMat("mutmp = ", mutmp, 1, ssdimITER);
//					Rprintf("dmvnorm = %f\n",dmvnorm(tUkx, mutmp, sig2mat, ssdimITER, log(det), scratch1, 0));
//					Rprintf("Ldmvnorm = %f\n",dmvnorm(tUkx, mutmp, sig2mat, ssdimITER, log(det), scratch1, 1));
//					RprintVecAsMat("nutmp", nutmp, 1, *C);
//					Rprintf("C = %d\n", *C);
					
				
//					Rprintf("dmultinom = %f\n", dmultinom(ytmp, nutmp, 1, *C, 0)); 

					mvnDens[h] = sbweight[h]*dmvnorm(tUkx, mutmp, sig2mat, ssdimITER, log(det), scratch1, 0);

//					Rprintf("mvnDens = %f\n", mvnDens[h]); 

					cwdensMix = cwdensMix + mvnDens[h];
					
					if(SiITER[j] == h+1){
						
						matrix_product(UmatITER, mutmp, Umu, *m, 1, ssdimITER);
						for(kk = 0; kk < *m; kk++) Umutheta[kk] = Umu[kk] + thetavec[kk];
							
//						loglike = log(mvnDens[h]) + dmultinom(ytmp, nutmp, 1, *C, 1);
						loglike = sbweight[h]*
								  dmvnorm(xtmp, Umutheta, sig2mat2, *m, log(det2), scratch1, 0)*
								  dmultinom(ytmp, nutmp, 1, *C, 0);

//						Rprintf("loglike = %f\n", log(loglike));
				
					}
				}
								
//				Rprintf("loglike = %f\n", loglike);
				for(h = 0; h < *N; h++){
				
					condwght[h] = mvnDens[h]/cwdensMix;

					fprintf(condwtsMCMC,"%10.9f\n", condwght[h]); 
					
				}

				fprintf(llikeMCMC,"%10.9f\n", loglike); 

			}
//			Rprintf("llike = %f\n", llike);
		}

*/


      //////////////////////////////////////////////////////////////////////////////////
      //																			  //
      // Save MCMC iterates															  //
      //																			  //
      //////////////////////////////////////////////////////////////////////////////////
      if((i > (*burn-1)) & ((i+1) % *thin == 0)){
//      	Rprintf("ii = %d\n", ii);
            
      	for(t = 0; t < *dim; t++){
      	  mu0[ii*(*dim) + t] = mu0_iter[t];
      	}
 
 
        for(d=0; d < ncp; d++){
          Sigma0[ii*ncp + d] = Sigma0_iter[d];
        }
 
        for(h = 0; h < *N; h++){
          sbw[ii*(*N) + h] = sbweight[h];
      	  for(t = 0; t < *dim; t++){
      	    mu[(ii*(*N) + h)*(*dim) + t] = mu_iter[h*(*dim) + t];
      	  }
      	  for(d = 0; d < ncp; d++){
      	    Sigma[(ii*(*N) + h)*(ncp) + d] = Sigma_iter[h*ncp + d];	
      	  }
      
      	}
      
 
      	for(j = 0; j < *nobs; j++){
      	  Si[ii*(*nobs) + j] = Si_iter[j];
      	  llike[ii*(*nobs) + j] = log(like_iter[j]);	
      	  fitted[ii*(*nobs) + j] = ispred[j];
      	}
      
      
      	nclus[ii] = nclus_iter;
      
      
      	ii = ii+1;
//      	Rprintf("ii = %d\n", ii);
      }
	}


	
	lpml_iter=0.0;
	for(j = 0; j < *nobs; j++){
//	  Rprintf("j = %d\n", j);
	  lpml_iter = lpml_iter - log(CPOinv[j]);
	}
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
