
data {
    int<lower=0> N;
    int<lower=0> P;
    array[N] real y;
    matrix[N,P] X;
    real<lower=0> sigma_alpha;
    real<lower=0> sigma_eps;
}
parameters{
    vector[P] beta;
    real<lower=0> tau;
    real alpha;
    real sigma;
}

model{
    // prior
    for (j in 1:P){
       beta[j] ~ double_exponential(0,tau);
       }
    tau ~ cauchy(0,1);
    alpha ~ normal(0,sigma_alpha);
    sigma ~ normal(0,sigma_eps);

    
    y ~ normal(alpha+X*beta,sigma);
      
} 

