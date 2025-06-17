data {
  int<lower=1> N;
  int<lower=1> N_miss;
  int<lower=1> N_obs;
  int<lower=1> D;
  vector[N_obs] y_obs;
  array[N] int<lower=0,upper=1> ry;
  matrix[N_obs,D] X_obs;
  matrix[N,D] X;
  real<lower=0> psi_1;
  real<lower=0> psi_2;
  real<lower=0> psi_3;
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
  vector[D] b_y1;
  vector[D] b_r;
  real<lower=0> sd1;
}
model {
  vector[N] mu_ry;
  vector[N_obs] mu_y_obs;

  b_y1 ~ normal(0, 1);
  b_r ~ normal(0, 1);
  sd1 ~ normal(0, 1);

  mu_ry = X * b_r;
  mu_y_obs = X_obs * b_y1;

  ry ~ bernoulli_logit(mu_ry);
  y_obs ~ normal(mu_y_obs, sd1);
}
generated quantities {
  vector[D] marg_b_y;
  real mu_y = 0;
  real y_ho;
  real sd0;
  {
    vector[D] b_y0;
    vector[N] mu_y_r1 = X * b_y1;
    vector[N] p_ry = inv_logit(X * b_r);
    b_y0[1] = psi_2 > 0 ? normal_rng(b_y1[1], psi_2^2 * abs(b_y1[1])) : b_y1[1];
    b_y0[2:D] = psi_1 > 0 ? b_y1[2:D] + normal_rng(0, psi_1) * b_y1[2:D] : b_y1[2:D];
    sd0 = psi_3 > 0 ? gamma_rng(sd1 / psi_3, 1 / psi_3) : sd1;
    vector[N] mu_y_r0 = X * b_y0;
    for (n in 1:N) {
      mu_y += ry[n] == 1 ? mu_y_r1[n] : mu_y_r0[n];
    }
    mu_y /= N;
    marg_b_y = b_y1 * mean(p_ry) + b_y0 * mean(1 - p_ry);
    real p_ho = inv_logit(X_ho * b_r);
    int r_ho = bernoulli_rng(p_ho);
    y_ho = r_ho ? normal_rng(X_ho * b_y1, sd1) : normal_rng(X_ho * b_y0, sd0);
  }
}
