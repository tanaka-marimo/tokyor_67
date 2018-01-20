functions {
  real ZINB_lpmf(int Y, real q, real mu, real phi) {
    if (Y == 0) {
      return log_sum_exp(
        bernoulli_lpmf(0 | q),
        bernoulli_lpmf(1 | q) + neg_binomial_2_log_lpmf(0 | mu, phi)
      );
    } else {
      return bernoulli_lpmf(1 | q) + neg_binomial_2_log_lpmf(Y | mu, phi);
    }
  }
}

data {
  int N;
  int Y[N];
  real X1[N];
  real X2[N];
  int S[N];
}

parameters {
  real a0;
  real a1;
  real a2;
  real c;
  real d0;
  real d1;
  real phi;
}

transformed parameters {
  real q[N];
  real mu[N];
  real c_new;
  real mu_new[N];

  for (n in 1:N)
    q[n] = inv_logit(d0 + d1*X1[n]);
  for (n in 1:N)
    mu[n] = a0 + a1*X1[n] + a2*X2[n] + c*S[n];
  c_new = 1 - (1/exp(c)); 
  
  for (n in 1:N)
    if(q[n]>=0.5){
      mu_new[n] = exp(mu[n]);
    } else {
      mu_new[n] = 0;
    }
}

model {
  for (n in 1:N)
    Y[n] ~ ZINB(q[n], mu[n], phi);
}
