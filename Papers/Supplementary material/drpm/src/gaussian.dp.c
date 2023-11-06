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

static void gibbs(double* y, int* n, int *N,
            double* m, double* v, double* A, double* A0, double *alpha, double *mh,
            int* niter, int *nburn, int *nthin,
            double* mu, double* sigma2, double* w, 
            int* z, double* theta, double* tau2,
            double* llike, double* fitted, double* lpml, double* waic ){


  int i, j, k, ii, nk;
  int nout = (*niter - *nburn)/(*nthin);

  // scratch vectors to update parameters
  double sumdens, cprob, uu, sumy=0.0, summu;
  double scr1[*N], dens[*N], sbw[*N], mnlike[*n], mnllike[*n], CPOinv[*n];
  double mstar, vstar, astar, bstar, ccount, cumccount, csbw;
  double llo, lln, llr, osig, nsig, otau, ntau;

  double _lpml, elppdWAIC, _theta, _tau2;
  double _mu[*N], _sigma2[*N], _w[*N], _like[*n], _fitted[*n];
  int _z[*n]; 
  
  /* Initialize variables */
  for(k = 0; k < *N; k++){ 
    _mu[k] = rnorm(0,1);
    _sigma2[k] = rgamma(1,1);
  }

  for(j = 0; j < *n; j++){
    _z[j] = rbinom(*N, 0.25);
    mnlike[j] = 0.0;
    mnllike[j] = 0.0;
    CPOinv[j] = 0.0;
  }
  
  _theta = 0.0; _tau2 = 1.0;

  // Since I update w first in the MCMC algorithm,
  // I don't need to initialize it
  
  Rprintf("Prior values: m0 = %.2f, s20 = %.2f\n \t A = %.2f, A0 = %.2f\n",
      	           *m, *v, *A, *A0);

  double candSIG = mh[0], candTAU = mh[1];

  
  ii = 0;

  for(i=0; i<*niter; i++){
    
    if(i % 5000 == 0) Rprintf("mcmc iter = %d\n", i);

    // 1) update stick breaking weights
    csbw = 1.0;
	for(k=0; k<*N; k++){
	  ccount = 0.0;
	  cumccount = 0.0;
	  for(j=0; j<*n; j++){

	    if(_z[j] == k+1){
		  ccount = ccount + 1.0;
	    }

		if(_z[j] > k+1){
	     cumccount = cumccount + 1.0;
		}
	  }

	  astar = 1.0 + ccount;
	  bstar = *alpha + cumccount;
					
	  sbw[k] = rbeta(astar,bstar);

	  if(k==(*N-1)) sbw[k] = 1.0;

	  if(k == 0){
		_w[k] = sbw[k];
	  }
	  if(k > 0){
	    csbw = csbw*(1-sbw[k-1]);
	    _w[k] = sbw[k]*csbw;
	  }			
	}
    
    // 2) update zi - component labels
    for(j=0; j<*n; j++){
      sumdens = 0.0;
      for(k=0; k<*N; k++){
        dens[k] = dnorm(y[j], _mu[k], sqrt(_sigma2[k]), 0)*_w[k]; 
        sumdens = sumdens + dens[k];
      }
//      RprintVecAsMat("dens", dens, 1, *N);
      for(k=0; k<*N; k++){
        scr1[k] = dens[k]/sumdens;
      }
//      RprintVecAsMat("scr1", scr1, 1, *N);

	  uu = runif(0.0,1.0);
	  cprob= 0.0;
      for(k = 0; k < *N; k++){
	    cprob = cprob + scr1[k];
		if (uu < cprob){
		  _z[j] = k+1;
		  break;
		}
	  }
    }
    
    // update mu
    summu = 0.0;
    for(k=0; k<*N; k++){
      sumy = 0.0;
      nk = 0;
      for(j=0; j<*n; j++){
        if(_z[j] == (k+1)){
          sumy = sumy + y[j];  
          nk = nk+1;
        }
      }

      vstar = 1/((nk/_sigma2[k]) + 1/(_tau2));
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
    for(k=0; k<*N; k++){
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
    vstar = 1.0/((*N/_tau2) + (1.0/(*v)));  
    mstar = vstar*((1/_tau2)*summu + (*m)/(*v));
    _theta = rnorm(mstar, sqrt(vstar));


    // update tau2
    otau = sqrt(_tau2);
    ntau = rnorm(otau, candTAU);
    if(ntau > 0.0){
      llo = 0.0;
      lln = 0.0;
       
      for(k=0; k<*N; k++){
        llo = llo + dnorm(mu[k], _theta, otau, 1);
        lln = lln + dnorm(mu[k], _theta, ntau, 1);
      }
        
      llo = llo + dunif(osig, 0, *A0, 1);
      lln = lln + dunif(nsig, 0, *A0, 1);
        
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
        _like[j] = 0.0;
        _fitted[j] = 0.0;
        for(k=0; k<*N; k++){
	      _like[j] = _like[j] + 
		                  _w[k]*
		                  dnorm(y[j], _mu[k], sqrt(_sigma2[k]), 0);
		  _fitted[j] = rnorm(_mu[_z[j]-1], sqrt(_sigma2[_z[j]-1]));
        }

        // These are needed for WAIC

        mnlike[j] = mnlike[j] + (_like[j])/((double) nout);
        mnllike[j] = mnllike[j] + log(_like[j])/((double) nout);

        CPOinv[j] = CPOinv[j] + (1/(double) nout)*(1/_like[j]);
      }
    
    
      // save  iterates
      theta[ii] = _theta;
      tau2[ii] = _tau2;

      for(k=0; k<*N; k++){
        mu[ii + nout*k] = _mu[k];
        sigma2[ii + nout*k] = _sigma2[k];
        w[ii + nout*k] = _w[k];
      }
      for(j=0; j<*n; j++){
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





SEXP GIBBS(SEXP y, SEXP n, SEXP N, SEXP m, SEXP v, SEXP A, SEXP A0, SEXP alpha,
	       SEXP mh, SEXP niter, SEXP nburn, SEXP nthin) {
  int nprot = 0;
  int _n = asInteger(n);
  int _N = asInteger(N);
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

  SEXP MU = PROTECT(allocMatrix(REALSXP, nout, _N)); nprot++;
  SEXP SIGMA2 = PROTECT(allocMatrix(REALSXP, nout, (_N))); nprot++; 
  SEXP W = PROTECT(allocMatrix(REALSXP, nout, (_N))); nprot++; 
  SEXP Z = PROTECT(allocMatrix(INTSXP, nout, _n)); nprot++; 
  SEXP LLIKE = PROTECT(allocMatrix(REALSXP, nout, _n)); nprot++;
  SEXP FITTED = PROTECT(allocMatrix(REALSXP, nout, _n)); nprot++;
  SEXP THETA = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP TAU2 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP WAIC = PROTECT(Rf_allocVector(REALSXP, 1)); nprot++;
  SEXP LPML = PROTECT(Rf_allocVector(REALSXP, 1)); nprot++;


  double *MUout, *SIGMA2out, *Wout, *THETAout, *TAU2out;
  double *LLIKEout, *WAICout, *LPMLout, *FITTEDout;
  int *Zout;

  MUout = REAL(MU);
  SIGMA2out = REAL(SIGMA2);
  Wout = REAL(W);
  Zout = INTEGER(Z);
  THETAout = REAL(THETA);
  TAU2out = REAL(TAU2);

  LLIKEout = REAL(LLIKE);
  FITTEDout = REAL(FITTED);
  WAICout = REAL(WAIC);
  LPMLout = REAL(LPML);

  GetRNGstate();

  gibbs(REAL(y), &_n, &_N, &_m, &_v, &_A, &_A0, &_alpha, REAL(mh), &_niter, &_nburn, &_nthin, 
        MUout, SIGMA2out, Wout, Zout, THETAout, TAU2out, LLIKEout, FITTEDout, LPMLout, WAICout);

  PutRNGstate();


  SEXP ans = PROTECT(allocVector(VECSXP, 10)); nprot++;
  SET_VECTOR_ELT(ans, 0, MU);
  SET_VECTOR_ELT(ans, 1, SIGMA2);
  SET_VECTOR_ELT(ans, 2, W);
  SET_VECTOR_ELT(ans, 3, Z);
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
  SET_STRING_ELT(nm, 2, mkChar("w"));
  SET_STRING_ELT(nm, 3, mkChar("Si"));
  SET_STRING_ELT(nm, 4, mkChar("theta"));
  SET_STRING_ELT(nm, 5, mkChar("tau2"));
  SET_STRING_ELT(nm, 6, mkChar("llike"));
  SET_STRING_ELT(nm, 7, mkChar("fitted"));
  SET_STRING_ELT(nm, 8, mkChar("lpml"));
  SET_STRING_ELT(nm, 9, mkChar("waic"));

  UNPROTECT(nprot);
  return(ans);
}



