library(parallel)
library(tidyverse)
library(here)

get_trial <- function(trial, rho_y_ry, mu_probit_ry) {

  # set.seed(trial)

  cover <- function(interval,true) {
    idx <- findInterval(true, c(-Inf,interval,Inf))
    if (idx == 2)
      return(1)
    else
      return(0)
  }

  gen_fit_summary <- function(fit, true_pars, model_nm) {
    par_sum <- fit$summary(c("marg_b_y", "y_ho"),mean, sd, ~quantile(.x, probs = c(0.025,0.975)))
    par_sum$true <- true_pars
    par_sum$bias <- par_sum$mean - par_sum$true
    par_sum$cover <- apply(par_sum[,c("2.5%","97.5%","true")],1, \(x) cover(x[1:2], x[3]))
    par_sum$model <- model_nm
    par_sum <- par_sum[,c("model", "variable", "bias", "cover", "sd")]
    return(par_sum)
  }

  ## rho_y_ry: Controls how MNAR y is
  ## mu_probit_ry: controls how missing y is
  gen_data <- function(n = 100, rho_y_ry = 0.5, mu_probit_ry = -0.5) {
    y_ry_Sigma <- matrix(c(1, rho_y_ry, rho_y_ry, 1),2,2)
    x2 <- rnorm(n, 0, 1)
    x1 <- rbinom(n, 1, 0.5)

    y_ry_noise <- MASS::mvrnorm(n, c(0,0), Sigma = y_ry_Sigma)
    y <- 1 + x1 + x2 + y_ry_noise[,1]
    ry_latent <- mu_probit_ry + x1 + x2 + y_ry_noise[,2]
    ry <- as.integer(ry_latent > 0)

    y_obs <- y[ry == 1]
    N_obs <- sum(ry)
    N_miss <- n - N_obs
    X <- cbind(1,x1,x2)
    X_obs <- X[ry == 1,]
    x1_ho <- rbinom(1, 1, 0.5)
    x2_ho <- rnorm(1)
    y_ry_noise_ho <- MASS::mvrnorm(1, c(0,0), Sigma = y_ry_Sigma)
    y_ho <- 1 + x1_ho + x2_ho + y_ry_noise_ho[1]
    stan_data <- list(
      N = n,
      N_obs = N_obs,
      N_miss = N_miss,
      D = ncol(X),
      y_obs = y_obs,
      ry = ry,
      X_obs = X_obs,
      X = X,
      mean_y = mean(y),
      y_ho = y_ho,
      X_ho = c(1, x1_ho, x2_ho)
    )
    return(stan_data)
  }

  heck <- cmdstan_model("inst/scripts/heckman.stan")
  pm <- cmdstan_model("inst/scripts/pattern-mix.stan")

  data <- gen_data(rho_y_ry = rho_y_ry)

  true_pars <- c(1,1,1,data$y_ho)
  heck_fit <- heck$sample(data = data, chains = 4)
  sum_heck <- gen_fit_summary(heck_fit, true_pars , "heckman")

  ## Assumes MNAR
  data_pm_p44 <- c(data, list(psi_1 = 0.4, psi_2 = 0.4, psi_3 = 0.4))
  pm_fit <- pm$sample(data = data_pm_p44, chains = 4)
  sum_p44 <- gen_fit_summary(pm_fit, true_pars , "pattern-mix-moderate-MNAR")

  ## Assumes more MNAR
  data_pm_11 <- c(data, list(psi_1 = 1, psi_2 = 1, psi_3 = 1))
  pm_fit <- pm$sample(data = data_pm_11, chains = 4)
  sum_11 <- gen_fit_summary(pm_fit, true_pars , "pattern-mix-large-MNAR")

  ## Assumes MAR/MCAR
  data_pm_00 <- c(data, list(psi_1 = 0, psi_2 = 0, psi_3 = 0))
  pm_fit <- pm$sample(data = data_pm_00, chains = 4)
  sum_mar <- gen_fit_summary(pm_fit, true_pars, "MAR/MCAR")

  total_res <- rbind(sum_heck, sum_p44, sum_11, sum_mar)
  total_res$trial <- trial
  total_res

}

nsim <- 2000
nc <- 50
cl <- makeCluster(nc)
clusterEvalQ(cl, {
  library(cmdstanr)
})
out <- parLapply(cl, 1:nsim, safely(get_trial), rho_y_ry = 0.1, mu_probit_ry = -0.5)
stopCluster(cl)

errors <- map_lgl(out, ~ !is.null(.x$error))
sum(errors)

out_success <- map(out, ~ .x$result)[!errors]

out_success <- out_success %>%
  bind_rows()

out_summary <- out_success %>%
  group_by(model, variable) %>%
  summarize(mbias = mean(bias), rmse = sqrt(mean(bias^2)), se = mean(sd), cover95 = mean(cover))
out_summary

write_csv(out_success, here("inst", "output", "simulation", "mnar_all_small.csv"))
write_csv(out_summary, here("inst", "output", "simulation", "mnar_summary_small.csv"))
