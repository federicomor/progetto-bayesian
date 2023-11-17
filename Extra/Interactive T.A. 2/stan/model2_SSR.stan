
data {
    int<lower=0> N;
    int<lower=0> P;
    array[N] real y;
    matrix[N,P] X;
    real<lower=0> sigma_alpha;
    real<lower=0> sigma_eps;
    real<lower=0> tau1 ;
    real<lower=0> tau2 ;

}
parameters{
    vector[P] beta;
    real alpha;
    real<lower=0> sigma;

}

model{
    // prior
    alpha ~ normal(0,sigma_alpha);
    sigma ~ normal(0,sigma_eps);
    for (j in 1:P) {
        target += log_mix(0.5,
            normal_lpdf(beta[j] | 0, tau1),
            normal_lpdf(beta[j] | 0, tau2)
        );
    }

    // Likelihood
    y ~ normal(X * beta + alpha, sigma);
} 


