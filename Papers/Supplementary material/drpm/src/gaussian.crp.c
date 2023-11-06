/****************************************************************************************
*
*
* My attempt to produce and MCMC algorith for a finite gaussian mixture using 
* the .Call function
*
*
****************************************************************************************/
#include "matrix.h"
#include "Rutil.h"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>

#include <math.h>
#include <stdio.h>
#include <stddef.h>

// inputs
// y - data vector
// n - length of y
// N - number of mixture components
// m - prior mean of the mean from mixture components
// v - prior variance of the mean from mixture components
// a - prior shape of the variance from mixture components
// b - prior rate of the variance from mixture components
// alpha - prior scale parameter of the DP 
// niter - number of MCMC iterates
// nburn - number of MCMC iterates to be discarded
// nthin - amount of thinning applied to the MCMC chain
//
// outputs
// mu - matrix containing MCMC draws for component means
// sigma2 - matrix containing MCMC draws for component variances
// pi - matrix containing MCMC draws for component weights
// z - matrix containing MCMC draws for component labels
// llike - matrix containing MCMC draws for log-likelihood values
// lpml - scalar containing the lpml value
// waic - scalar containing the waic value

static void gibbs_crp(double* y, int* n, 
            double* m, double* v, double* A, double* A0, double* alpha, double* mh,
            int* niter, int *nburn, int *nthin,
            double* mu, double* sigma2, 
            int* z, double* theta, double* tau2, int* nclus,
            double* llike, double* fitted, double* lpml, double* waic ){


  int i, j, jj, k, ii;
  int nout = (*niter - *nburn)/(*nthin);
  double M = *alpha;

  // scratch vectors to update parameters
  double denph, cprob, uu, sumy=0.0, summu;
  double ph[*n], mnlike[*n], mnllike[*n], CPOinv[*n];
  double mstar, vstar, mudraw, sigdraw, maxph;
  double llo, lln, llr, osig, nsig, otau, ntau, auxs;

  double _lpml, elppdWAIC, _theta, _tau2;
  double _mu[*n], _sigma2[*n], _like[*n], _fitted[*n];
  int _z[*n], nh[*n], _nclus=0, iaux; 
  
  /* Initialize variables */
  for(j = 0; j < *n; j++){
    _z[j] = rbinom(1, 0.25)+1;
    _mu[j] = 0;
    _sigma2[j] = 1;
    mnlike[j] = 0.0;
    mnllike[j] = 0.0;
    CPOinv[j] = 0.0;
    nh[j] = 0;
  }
  
  // Initial enumeration of number of observations per cluster;
  for(j = 0; j < *n; j++){
    for(jj = 0; jj < *n; jj++){
	  if(_z[j] == jj+1) nh[jj] = nh[jj] + 1;
	}
  }
  // Initialize the number of clusters	
  for(j = 0; j < *n; j++){
    if(nh[j] > 0) _nclus = _nclus + 1;
  }

  _theta = 0.0; _tau2 = 1.0;



  // Since I update w first in the MCMC algorithm,
  // I don't need to initialize it
  
  Rprintf("Prior values: m0 = %.2f, s20 = %.2f\n \t A = %.2f, A0 = %.2f\n",
      	           *m, *v, *A, *A0);

  double candSIG = mh[0], candTAU = mh[1];

  
  ii = 0;

  for(i=0; i<*niter; i++){
    
    if(i % 10000 == 0) Rprintf("mcmc iter = %d\n", i);


    // Neal's algorithm 8 that integrates out the random probability measure
    for(j=0; j<*n; j++){
      if(nh[_z[j]-1] > 1){
	    // Observation belongs to a non-singleton ...
        nh[_z[j]-1] = nh[_z[j]-1] - 1;
	  }else{
	    iaux = _z[j];
//		Rprintf("iaux = %d\n", iaux);
		if(iaux < _nclus){
				
		  // Need to relabel clusters. I will do this by swapping cluster labels
		  // _z[j] and _nclus along with cluster specific parameters;
				
		  // All members of last cluster will be assigned subject j's label
		  for(jj=0; jj<*n; jj++){
		    if(_z[jj] == _nclus){
			  _z[jj] = iaux; 
			}
		  }
		  _z[j] = _nclus;

		  // The following steps swaps order of cluster specific parameters
		  // so that the newly labeled subjects from previous step retain
		  // their correct cluster specific parameters
		  auxs = _sigma2[iaux-1];
		  _sigma2[iaux-1] = _sigma2[_nclus-1];
		  _sigma2[_nclus-1] = auxs;

		  auxs = _mu[iaux-1];
		  _mu[iaux-1] = _mu[_nclus-1];
		  _mu[_nclus-1] = auxs;

		  // the number of members in cluster is also swapped with the last
		  nh[iaux-1] = nh[_nclus-1];
		  nh[_nclus-1] = 1;
			
		}
			
		// Now remove the ith obs and last cluster;
		nh[_nclus-1] = nh[_nclus-1] - 1;
		_nclus = _nclus - 1;
			
	  }
	  
	  ///////////////////////////////////
	  //	
	  // Begin the cluster probabilities
	  //
	  //////////////////////////////////
	  for(k=0; k<_nclus; k++){

	    ph[k] = dnorm(y[j], _mu[k], sqrt(_sigma2[k]), 1) +
				  log((double) nh[k]) - 
				  log(*n - 1 + M);
		
	  }

      // Need to consider a singleton
      mudraw = rnorm(_theta, sqrt(_tau2));
      sigdraw = runif(0, *A);
			
	  ph[_nclus] = dnorm(y[j], mudraw, sigdraw, 1) +
			         log(M) -  
			         log(*n - 1 + M);


	  maxph = ph[0];
	  for(k=1; k<_nclus+1; k++){
	    if(ph[k] > maxph) maxph = ph[k];
	  }
	  
	  denph = 0.0;
	  for(k = 0; k < _nclus+1; k++){
	    ph[k] = exp(ph[k] - maxph);
		denph = denph + ph[k];
	  }

	  for(k = 0; k < _nclus+1; k++){
	    ph[k] = ph[k]/denph;
	  }
	
	  uu = runif(0.0,1.0);
		
	  cprob= 0.0;;
	  for(k = 0; k < _nclus+1; k++){
	    cprob = cprob + ph[k];
		if (uu < cprob){
		  iaux = k+1;
		  break;	
		}
	  }		

	  if(iaux <= _nclus){
	    _z[j] = iaux;
		nh[_z[j]-1] = nh[_z[j]-1] + 1;
	  }else{
	    _nclus = _nclus + 1;
		_z[j] = _nclus;
		nh[_z[j]-1] = 1;			
				
		_mu[_z[j]-1] = mudraw;
		_sigma2[_z[j]-1] = sigdraw*sigdraw;

	  }
	}


    
    // update mu
    summu = 0.0;
    for(k=0; k<_nclus; k++){
      sumy = 0.0;
      for(j=0; j<*n; j++){
        if(_z[j] == (k+1)){
          sumy = sumy + y[j];  
        }
      }

      vstar = 1/((nh[k]/_sigma2[k]) + 1/(_tau2));
      mstar = vstar*((1.0/_sigma2[k])*sumy + (_theta)/(_tau2));
      
      _mu[k] = rnorm(mstar, sqrt(vstar));
      summu = summu + _mu[k];
    }


    // update sigma2
/*    for(k=0; k<*N; k++){
      ssq=0.0, nk=0;
      for(j=0; j<*n; j++){
        if(_z[j] == (k+1)){
          ssq = ssq + (y[j]-_mu[k])*(y[j]-_mu[k]);
          nk = nk + 1;
        }
      }
      astar = 0.5*(nk) + *a;
      bstar = 0.5*ssq + *b;
      _sigma2[k] = 1/rgamma(astar, 1/bstar);
    }
*/

    // update sigma2
    for(k=0; k<_nclus; k++){
      osig = sqrt(_sigma2[k]);
      nsig = rnorm(osig, candSIG);
      if(nsig > 0){
        llo = 0.0;
        lln = 0.0;
       
        for(j=0; j<*n; j++){
          if(_z[j] == (k+1)){
            llo = llo + dnorm(y[j], mu[k], osig, 1);
            lln = lln + dnorm(y[j], mu[k], nsig, 1);
          }
        }
        llo = llo + dunif(osig, 0, *A, 1);
        lln = lln + dunif(nsig, 0, *A, 1);
        
        llr = lln - llo;
        uu = runif(0,1);
        if(llr > log(uu)){
          _sigma2[k] = nsig*nsig;  
        }
      }
    }


    // update theta
    vstar = 1.0/((_nclus/_tau2) + (1.0/(*v)));  
    mstar = vstar*((1/_tau2)*summu + (*m)/(*v));
    _theta = rnorm(mstar, sqrt(vstar));


    // update tau2
    otau = sqrt(_tau2);
    ntau = rnorm(otau, candTAU);

    if(ntau > 0.0){
      llo = 0.0;
      lln = 0.0;
       
      for(k=0; k<_nclus; k++){
        llo = llo + dnorm(_mu[k], _theta, otau, 1);
        lln = lln + dnorm(_mu[k], _theta, ntau, 1);
      }

      llo = llo + dunif(otau, 0, *A0, 1);
      lln = lln + dunif(ntau, 0, *A0, 1);

      llr = lln - llo;
      uu = runif(0,1);
      if(llr > log(uu)){
        _tau2 = ntau*ntau;  
      }
    }
    
    // evaluate likelihood and keep iterates and produce fitted values
    if((i > (*nburn-1)) & ((i+1) % *nthin ==0)){
    
      // evaluate likelihood to Compute the CPO and lpml
      for(j=0; j<*n; j++){

        _like[j] = dnorm(y[j], _mu[_z[j]-1], sqrt(_sigma2[_z[j]-1]), 0);
		_fitted[j] = rnorm(_mu[_z[j]-1], sqrt(_sigma2[_z[j]-1]));


        // These are needed for WAIC

        mnlike[j] = mnlike[j] + (_like[j])/((double) nout);
        mnllike[j] = mnllike[j] + log(_like[j])/((double) nout);

        CPOinv[j] = CPOinv[j] + (1/(double) nout)*(1/_like[j]);
      }
    
    
      // save  iterates
      theta[ii] = _theta;
      tau2[ii] = _tau2;
      nclus[ii] = _nclus;
      
      for(j=0; j<*n; j++){
        mu[ii + nout*j] = _mu[_z[j]-1];
        sigma2[ii + nout*j] = _sigma2[_z[j]-1];
        z[ii + nout*j] = _z[j];
        llike[ii + nout*j] = log(_like[j]);
        fitted[ii + nout*j] = _fitted[j];
      }

      ii = ii + 1;
    }
  }

  _lpml=0.0;
  for(j=0; j<*n; j++){
    _lpml = _lpml + log(1/CPOinv[j]);
  }
  lpml[0] = _lpml;
  
  ////////////////////////////////////////////////////////////////////////////////////////////
  // Computing WAIC  (see Gelman article in lit review folder)
  ////////////////////////////////////////////////////////////////////////////////////////////
  elppdWAIC = 0.0;
  for(j=0; j<*n; j++){
    elppdWAIC = elppdWAIC + (2*mnllike[j] - log(mnlike[j]));
  }
  waic[0] = -2*elppdWAIC;
  
}



SEXP GIBBS_CRP(SEXP y, SEXP n, SEXP m, SEXP v, SEXP A, SEXP A0, SEXP alpha,
	       SEXP mh, SEXP niter, SEXP nburn, SEXP nthin) {
  int nprot = 0;
  int _n = asInteger(n);
  int _niter = asInteger(niter);
  int _nburn = asInteger(nburn);
  int _nthin = asInteger(nthin);
  double _m = asReal(m);
  double _v = asReal(v);
  double _A = asReal(A);
  double _A0 = asReal(A0);
  double _alpha = asReal(alpha); 

  double nout = (_niter-_nburn)/_nthin;
  
  y = PROTECT(coerceVector(y, REALSXP)); nprot++;
  mh = PROTECT(coerceVector(mh, REALSXP)); nprot++;

  SEXP MU = PROTECT(allocMatrix(REALSXP, nout, _n)); nprot++;
  SEXP SIGMA2 = PROTECT(allocMatrix(REALSXP, nout, _n)); nprot++; 
  SEXP Z = PROTECT(allocMatrix(INTSXP, nout, _n)); nprot++; 
  SEXP NCLUS = PROTECT(allocMatrix(INTSXP, nout, 1)); nprot++; 
  SEXP LLIKE = PROTECT(allocMatrix(REALSXP, nout, _n)); nprot++;
  SEXP FITTED = PROTECT(allocMatrix(REALSXP, nout, _n)); nprot++;
  SEXP THETA = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP TAU2 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP WAIC = PROTECT(Rf_allocVector(REALSXP, 1)); nprot++;
  SEXP LPML = PROTECT(Rf_allocVector(REALSXP, 1)); nprot++;


  double *MUout, *SIGMA2out, *THETAout, *TAU2out;
  double *LLIKEout, *WAICout, *LPMLout, *FITTEDout;
  int *Zout, *NCLUSout;

  MUout = REAL(MU);
  SIGMA2out = REAL(SIGMA2);
  Zout = INTEGER(Z);
  NCLUSout = INTEGER(NCLUS);
  THETAout = REAL(THETA);
  TAU2out = REAL(TAU2);

  LLIKEout = REAL(LLIKE);
  FITTEDout = REAL(FITTED);
  WAICout = REAL(WAIC);
  LPMLout = REAL(LPML);

  GetRNGstate();

  gibbs_crp(REAL(y), &_n, &_m, &_v, &_A, &_A0, &_alpha, REAL(mh), &_niter, &_nburn, &_nthin, 
        MUout, SIGMA2out, Zout, THETAout, TAU2out, NCLUSout, LLIKEout, FITTEDout, LPMLout, WAICout);

  PutRNGstate();


  SEXP ans = PROTECT(allocVector(VECSXP, 10)); nprot++;
  SET_VECTOR_ELT(ans, 0, MU);
  SET_VECTOR_ELT(ans, 1, SIGMA2);
  SET_VECTOR_ELT(ans, 2, Z);
  SET_VECTOR_ELT(ans, 3, NCLUS);
  SET_VECTOR_ELT(ans, 4, THETA);
  SET_VECTOR_ELT(ans, 5, TAU2);
  SET_VECTOR_ELT(ans, 6, LLIKE);
  SET_VECTOR_ELT(ans, 7, FITTED);
  SET_VECTOR_ELT(ans, 8, LPML);
  SET_VECTOR_ELT(ans, 9, WAIC);


  SEXP nm = allocVector(STRSXP, 10);
  setAttrib(ans, R_NamesSymbol, nm);
  SET_STRING_ELT(nm, 0, mkChar("mu"));
  SET_STRING_ELT(nm, 1, mkChar("sigma2"));
  SET_STRING_ELT(nm, 2, mkChar("Si"));
  SET_STRING_ELT(nm, 3, mkChar("nclus"));
  SET_STRING_ELT(nm, 4, mkChar("theta"));
  SET_STRING_ELT(nm, 5, mkChar("tau2"));
  SET_STRING_ELT(nm, 6, mkChar("llike"));
  SET_STRING_ELT(nm, 7, mkChar("fitted"));
  SET_STRING_ELT(nm, 8, mkChar("lpml"));
  SET_STRING_ELT(nm, 9, mkChar("waic"));

  UNPROTECT(nprot);
  return(ans);
}



