
    data {
        int<lower=0> dim;
        matrix[dim, dim] cov_chol;
    }
    
    parameters {
        vector[dim] x;
    }
    
    model {
        vector[dim] mu = rep_vector(0, dim);
        x ~ multi_normal_cholesky(mu, cov_chol);
    }

