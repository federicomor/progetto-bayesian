
data {
    
    int<lower=0> K; // number of variables
    int<lower=0> J; // number of groups
    int<lower=0> N; // number of observations

    real<lower=0> tau1; // hyperparameter for the spike
    real<lower=0> tau2; // hyperparameter for the slab
    real<lower=0> a; // hyperparameter for Beta distribution
    real<lower=0> b; // hyperparameter for Beta distribution
}

parameters { }

transformed parameters { }

model {
    // Priors
    for (j in 1:J) {
        p0[j] ~ beta(a, b); // prior for p0
        for (k in 1:K) {
            target += log_mix(p0[j],
                              normal_lpdf(p0[j]| 0, tau1), 
                              normal_lpdf(p0[j] | 0, tau2)); 
        }
    }

    // Likelihood

}



