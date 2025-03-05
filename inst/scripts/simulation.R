library(tidyverse)
library(mice)
library(brms)
library(parallel)
library(here)


get_simulation <- function(trial) {

  get_bayesian_results <- function(trial, dat_train, dat_test, beta) {
    set.seed(trial)
    # the missing dat_traina for x2 should have a prior and be distributed as
    # gaussian with an intercept, then it is used in model
    # this is very different from FCS where each variable needs to be included in
    # the other's imputation model
    bform <- bf(y | mi() ~ x1 + mi(x2)) +
      bf(x2 | mi() ~ 1) + set_rescor(FALSE)
    full_bayes <- brm(bform, data = dat_train, cores = 1)
    full_bayes_fixef <- fixef(full_bayes)
    # remove x2 intercept from bayesian model
    full_bayes_fixef <- as.data.frame(full_bayes_fixef[-2, ])
    tidy_fits <- data.frame(
      term = c("(Intercept)", "x1B", "x2"),
      estimate = full_bayes_fixef$Estimate,
      std.error = full_bayes_fixef$Est.Error,
      lower = full_bayes_fixef$Q2.5,
      upper = full_bayes_fixef$Q97.5,
      true = beta,
      trial = trial
    ) %>%
      mutate(cover = lower < true & true < upper) %>%
      mutate(error = true - estimate) %>%
      mutate(method = "bayesian")

    draws <- as_draws_df(full_bayes) %>%
      select(contains("Ymi_x2"))
    x2_miss_index <- is.na(dat_train$x2)
    x2_true_miss <- dat_train$x2_true[x2_miss_index]
    x2_errors <- sapply(seq_len(length(x2_true_miss)), function(x) {
      draws[[x]] - x2_true_miss[x]
    })
    x2_bias <- mean(x2_errors)
    x2_mspe <- mean(x2_errors^2)

    draws <- as_draws_df(full_bayes) %>%
      select(contains("Ymi_y"))
    y_miss_index <- is.na(dat_train$y)
    y_true_miss <- dat_train$y_true[y_miss_index]
    y_errors <- sapply(seq_len(length(y_true_miss)), function(x) {
      draws[[x]] - y_true_miss[x]
    })
    y_bias <- mean(y_errors)
    y_mspe <- mean(y_errors^2)

    imp_stat <- tibble(x2_bias = x2_bias, y_bias = y_bias,
                       x2_mspe = x2_mspe, y_mspe = y_mspe,
                       method = "bayesian", trial = trial)

    pred_out <- predict(full_bayes, dat_test)[, , "y"]
    preds_fits <- data.frame(
      estimate = pred_out[["Estimate"]],
      std.error = pred_out[["Est.Error"]],
      lower = pred_out[["Q2.5"]],
      upper = pred_out[["Q97.5"]],
      true = dat_test$y,
      trial = trial
    ) %>%
      mutate(cover = lower < true & true < upper) %>%
      mutate(error = true - estimate) %>%
      mutate(method = "bayesian")

    list(fixed = tidy_fits, imp_stat = imp_stat, pred = preds_fits)
  }

  get_lm_results <- function(trial, dat_train, dat_test, beta) {
    lmod <- lm(y ~ x1 + x2, dat_train)
    tidy_fits <- tidy(lmod, conf.int = TRUE, conf.level = 0.95) %>%
      select(term, estimate, std.error, lower = conf.low, upper = conf.high) %>%
      mutate(true = beta, trial = trial) %>%
      mutate(cover = lower < true & true < upper) %>%
      mutate(error = true - estimate) %>%
      mutate(method = "cc")

    pred_out <- predict(lmod, newdata = dat_test, interval = "prediction", level = 0.95)
    preds_fits <- data.frame(
      estimate = pred_out[, "fit"],
      std.error = (pred_out[, "upr"] - pred_out[, "lwr"]) / (2 * qt(0.975, df = n - length(beta))),
      lower = pred_out[, "lwr"],
      upper = pred_out[, "upr"]
    )
    preds_fits <- preds_fits %>%
      mutate(true = dat_test$y, trial = trial) %>%
      mutate(cover = lower < true & true < upper) %>%
      mutate(error = true - estimate) %>%
      mutate(method = "cc")


    imp_stat <- tibble(x2_bias = NA, y_bias = NA,
                       x2_mspe = NA, y_mspe = NA,
                       method = "cc", trial = trial)

    list(fixed = tidy_fits, imp_stat = imp_stat, pred = preds_fits)
  }

  get_imputed_results <- function(trial, method, multiple = TRUE, formulas, dat_train,
                                  dat_test, beta, donors = 3L, matchtype = 1L, maxit = 5) {

    set.seed(trial)

    if (multiple) {
      m <- 100
    } else {
      m <- 1
    }

    # if (method == "norm.boot") {
    #   browser()
    # }
    imps <- mice(
      dat_train,
      m = m,
      formulas = formulas,
      method = method,
      printFlag = FALSE,
      donors = donors,
      matchtype = matchtype,
      maxit = maxit
    )

    x2_miss_index <- is.na(imps$data$x2)
    x2_errors <- sapply(m, function(x) {
      dat_train$x2_true[x2_miss_index] - complete(imps, m)$x2[x2_miss_index]
      })
    x2_bias <- mean(x2_errors)
    x2_mspe <- mean(x2_errors^2)

    y_miss_index <- is.na(imps$data$y)
    y_errors <- sapply(m, function(x) {
      dat_train$y_true[y_miss_index] - complete(imps, m)$y[y_miss_index]
    })
    y_bias <- mean(y_errors)
    y_mspe <- mean(y_errors^2)

    if (method == "pmm") {
      method <- paste0(method, "_", donors, "D", "_", "type", matchtype)
    }

    if (multiple) {
      method <- paste0(method, "_multiple")
    }

    lmods <- with(imps, lm(y ~ x1 + x2))
    fits <- pool(lmods)
    tidy_fits <- tidy(fits, conf.int = TRUE, conf.level = 0.95) %>%
      select(term, estimate, std.error, lower = conf.low, upper = conf.high) %>%
      mutate(true = beta, trial = trial) %>%
      mutate(cover = lower < true & true < upper) %>%
      mutate(error = true - estimate) %>%
      mutate(method = method)

    imp_stat <- tibble(x2_bias = x2_bias, y_bias = y_bias,
                       x2_mspe = x2_mspe, y_mspe = y_mspe,
                       method = method, trial = trial)

    if (multiple) {
      preds <- lapply(lmods$analyses, function(x) {
        pred_fits <- predict(x, newdata = dat_test, se.fit = TRUE)
        data.frame(fit = pred_fits$fit, var = pred_fits$se.fit^2 + pred_fits$residual.scale^2)
      })
      Q <- map_dbl(preds, ~ .x$fit)
      U <- map_dbl(preds, ~ .x$var)
      n <- NROW(dat_train)
      pred_out <- pool.scalar(Q, U, n)
      df_out <- pred_out$df
      se_out <- sqrt(pred_out$t)
      crit <- qt(0.975, df = df_out)
      preds_fits <- data.frame(
        estimate = pred_out$qbar,
        std.error = se_out,
        lower = pred_out$qbar - crit * se_out,
        upper = pred_out$qbar + crit * se_out
      )
    } else {
      pred_out <- predict(lmods$analyses[[1]], newdata = dat_test, interval = "prediction", level = 0.95)
      preds_fits <- data.frame(
        estimate = pred_out[, "fit"],
        std.error = (pred_out[, "upr"] - pred_out[, "lwr"]) / (2 * qt(0.975, df = n - length(beta))),
        lower = pred_out[, "lwr"],
        upper = pred_out[, "upr"]
      )
    }
    preds_fits <- preds_fits %>%
      mutate(true = dat_test$y, trial = trial) %>%
      mutate(cover = lower < true & true < upper) %>%
      mutate(error = true - estimate) %>%
      mutate(method = method)
    list(fixed = tidy_fits, imp_stat = imp_stat, pred = preds_fits)
  }

  # functions first

  set.seed(trial)
  # specify data
  n_obs <- 99
  n_pred <- 1
  n <- n_obs + n_pred
  x1 <- rep(c("A", "B"), each = n / 2)
  x2 <- rnorm(n, mean = 0)
  beta <- c(1, 1, 1)
  y <- beta[1] + beta[2] * (x1 == "B") + beta[3] * x2 + rnorm(n, sd = 1)
  dat <- data.frame(x1, x2, y)
  dat <- dat[sample(1:n), , drop = FALSE]
  dat_train <- dat[1:n_obs, , drop = FALSE]
  dat_test <- dat[(n_obs + 1):n, , drop = FALSE]

  get_binom <- function(x1, p1 = 0.1, p2 = 0.8) {
    if (x1 == "A") {
      val <- rbinom(1, 1, p1)
    } else {
      val <- rbinom(1, 1, p2)
    }
    val
  }

  dat_train <- dat_train %>%
    mutate(x2_true = x2) %>%
    mutate(y_true = y) %>%
    mutate(x2_BINOM = map_dbl(x1, get_binom, p1 = 0.35, p2 = 0.70)) %>%
    mutate(y_BINOM = map_dbl(x1, get_binom, p1 = 0.70, p2 = 0.35)) %>%
    mutate(x2 = if_else(x2_BINOM == 1, NA, x2)) %>%
    mutate(y = if_else(y_BINOM == 1, NA, y))

  formulas <- list(
    y ~ x1 + x2,
    x2 ~ x1 + y
  )

  formulas <- name.formulas(formulas)


  imp_mean <- get_imputed_results(trial, "mean", multiple = FALSE, formulas, dat_train, dat_test, beta = beta)
  imp_reg <- get_imputed_results(trial, "norm.predict", multiple = FALSE, formulas, dat_train, dat_test, beta = beta)
  imp_streg <- get_imputed_results(trial, "norm.nob", multiple = FALSE, formulas, dat_train, dat_test, beta = beta)
  imp_streg_mult <- get_imputed_results(trial, "norm.nob", multiple = TRUE, formulas, dat_train, dat_test, beta = beta)
  imp_boot <- get_imputed_results(trial, "norm.boot", multiple = FALSE, formulas, dat_train, dat_test, beta = beta)
  imp_boot_mult <- get_imputed_results(trial, "norm.boot", multiple = TRUE, formulas, dat_train, dat_test, beta = beta)
  imp_norm <- get_imputed_results(trial, "norm", multiple = FALSE, formulas, dat_train, dat_test, beta = beta)
  imp_norm_mult <- get_imputed_results(trial, "norm", multiple = TRUE, formulas, dat_train, dat_test, beta = beta)
  imp_lasso <- get_imputed_results(trial, "lasso.norm", multiple = FALSE, formulas, dat_train, dat_test, beta = beta)
  imp_lasso_mult <- get_imputed_results(trial, "lasso.norm", multiple = TRUE, formulas, dat_train, dat_test, beta = beta)
  imp_rf <- get_imputed_results(trial, "rf", multiple = FALSE, formulas, dat_train, dat_test, beta = beta)
  imp_rf_mult <- get_imputed_results(trial, "rf", multiple = TRUE, formulas, dat_train, dat_test, beta = beta)
  imp_pmm_1D_type0 <- get_imputed_results(trial, "pmm", multiple = FALSE, formulas, dat_train, dat_test, beta = beta, donors = 1L, matchtype = 0)
  imp_pmm_3D_type0 <- get_imputed_results(trial, "pmm", multiple = FALSE, formulas, dat_train, dat_test, beta = beta, matchtype = 0)
  imp_pmm_mult_3D_type0 <- get_imputed_results(trial, "pmm", multiple = TRUE, formulas, dat_train, dat_test, beta = beta, matchtype = 0)
  imp_pmm_1D <- get_imputed_results(trial, "pmm", multiple = FALSE, formulas, dat_train, dat_test, beta = beta, donors = 1L)
  imp_pmm_mult_1D <- get_imputed_results(trial, "pmm", multiple = TRUE, formulas, dat_train, dat_test, beta = beta, donors = 1L)
  imp_pmm_3D <- get_imputed_results(trial, "pmm", multiple = FALSE, formulas, dat_train, dat_test, beta = beta)
  imp_pmm_mult_3D <- get_imputed_results(trial, "pmm", multiple = TRUE, formulas, dat_train, dat_test, beta = beta)
  bayes <- get_bayesian_results(trial, dat_train, dat_test, beta = beta)
  cc <- get_lm_results(trial, dat_train, dat_test, beta = beta)

  fixed <- bind_rows(imp_mean$fixed, imp_reg$fixed, imp_streg$fixed, imp_streg_mult$fixed,
                     imp_boot$fixed, imp_boot_mult$fixed,
                     imp_norm$fixed, imp_norm_mult$fixed, imp_lasso$fixed, imp_lasso_mult$fixed,
                     imp_rf$fixed, imp_rf_mult$fixed,
                     imp_pmm_1D_type0$fixed, imp_pmm_3D_type0$fixed, imp_pmm_mult_3D_type0$fixed,
                     imp_pmm_1D$fixed, imp_pmm_mult_1D$fixed, imp_pmm_3D$fixed, imp_pmm_mult_3D$fixed,
                     bayes$fixed, cc$fixed)

  imp_stat <- bind_rows(imp_mean$imp_stat, imp_reg$imp_stat, imp_streg$imp_stat, imp_streg_mult$imp_stat,
                     imp_boot$imp_stat, imp_boot_mult$imp_stat,
                     imp_norm$imp_stat, imp_norm_mult$imp_stat, imp_lasso$imp_stat, imp_lasso_mult$imp_stat,
                     imp_rf$imp_stat, imp_rf_mult$imp_stat,
                     imp_pmm_1D_type0$imp_stat, imp_pmm_3D_type0$imp_stat, imp_pmm_mult_3D_type0$imp_stat,
                     imp_pmm_1D$imp_stat, imp_pmm_mult_1D$imp_stat, imp_pmm_3D$imp_stat, imp_pmm_mult_3D$imp_stat,
                     bayes$imp_stat, cc$imp_stat)

  pred <- bind_rows(imp_mean$pred, imp_reg$pred, imp_streg$pred, imp_streg_mult$pred,
                    imp_boot$pred, imp_boot_mult$pred,
                    imp_norm$pred, imp_norm_mult$pred, imp_lasso$pred, imp_lasso_mult$pred,
                    imp_rf$pred, imp_rf_mult$pred,
                    imp_pmm_1D_type0$pred, imp_pmm_3D_type0$pred, imp_pmm_mult_3D_type0$pred,
                    imp_pmm_1D$pred, imp_pmm_mult_1D$pred, imp_pmm_3D$pred, imp_pmm_mult_3D$pred,
                    bayes$pred, cc$pred)

  list(fixed = fixed, imp_stat = imp_stat, pred = pred)

}

nsim <- 2004
nc <- 24
cl <- makeCluster(nc)
clusterEvalQ(cl, {
  library(tidyverse)
  library(mice)
  library(brms)
})
out <- parLapply(cl, 1:nsim, safely(get_simulation))
stopCluster(cl)

# number of errors
errors <- map_lgl(out, ~ !is.null(.x$error))
sum(errors)

out_success <- map(out, ~ .x$result)[!errors]

fixed <- bind_rows(map(out_success, ~ .x$fixed))
fixed_summary <- fixed %>%
  group_by(term, method) %>%
  summarize(mbias = mean(error), mse = mean(error^2), mvar = mean(std.error^2), cover = mean(cover), .groups = "drop")

imp_stat <- bind_rows(map(out_success, ~ .x$imp_stat))
imp_stat_summary <- imp_stat %>%
  group_by(method) %>%
  summarize(x2_bias = mean(x2_bias), x2_mspe = mean(x2_mspe),
            y_bias = mean(y_bias), y_mspe = mean(y_mspe))

pred <- bind_rows(map(out_success, ~ .x$pred))
pred_summary <- pred %>%
  group_by(method) %>%
  summarize(mbias = mean(error), mspe = mean(error^2), mvar = mean(std.error^2), cover = mean(cover))

write_csv(fixed, here("inst", "output", "simulation", "fixed_all.csv"))
write_csv(imp_stat, here("inst", "output", "simulation", "imp_stat.csv"))
write_csv(pred, here("inst", "output", "simulation", "pred_all.csv"))
write_csv(fixed_summary, here("inst", "output", "simulation", "fixed_summary.csv"))
write_csv(imp_stat_summary, here("inst", "output", "simulation", "imp_stat_summary.csv"))
write_csv(pred_summary, here("inst", "output", "simulation", "pred_summary.csv"))

