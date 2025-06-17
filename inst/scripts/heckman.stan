data {
  int<lower=1> N;
  int<lower=1> N_miss;
  int<lower=1> N_obs;
  int<lower=1> D;
  vector[N_obs] y_obs;
  array[N] int<lower=0,upper=1> ry;
  matrix[N_obs,D] X_obs;
  matrix[N,D] X;
  row_vector[D] X_ho;
}
transformed data {
  array[N_obs] int n_obs;
  array[N_miss] int n_miss;
  {
    int i = 1;
    int j = 1;
    for (n in 1:N) {
      if (ry[n] == 1) {
        n_obs[i] = n;
        i = i + 1;
      } else {
        n_miss[j] = n;
        j = j + 1;
      }
    }
  }
}
parameters {
  real<lower=-1,upper=1> rho;
  vector[D] b_y;
  vector[D] b_r;
  real<lower=0> sd1;
}
model {
  vector[N] mu_ry;
  vector[N_obs] mu_y_obs;

  b_y ~ normal(0, 1);
  b_r ~ normal(0, 1);
  sd1 ~ normal(0, 1);

  mu_ry = X * b_r;
  mu_y_obs = X_obs * b_y;

  for (n in 1:N_miss) {
    target += (log(Phi(-mu_ry[n_miss[n]])));
  }
  for (n in 1:N_obs) {
    target += log(Phi(sqrt(1 - rho^2)^(-1)*(mu_ry[n_obs[n]]
                               + rho / sd1
                               * (y_obs[n] - mu_y_obs[n]))));
  }
  y_obs ~ normal(mu_y_obs, sd1);
}
generated quantities {
  real mu_y = mean(X * b_y);
  real y_ho = normal_rng(X_ho * b_y, sd1);
  vector[D] marg_b_y = b_y;
}
