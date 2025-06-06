library(tidyverse)
library(mice)
library(parallel)
library(here)

get_trial <- function(trial) {


  get_single_trial <- function(trial, include_y, mice_method) {

    set.seed(trial)

    if (include_y) {
      part3 <- "Yes"
    } else {
      part3 <- "No"
    }
    method <- str_c(mice_method, part3, sep = "_")
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


    if (include_y) {
      formulas <- list(
        y ~ x1 + x2,
        x2 ~ x1 + y
      )
    } else {
      formulas <- list(
        y ~ x1 + x2,
        x2 ~ x1
      )
    }
    formulas <- name.formulas(formulas)


    if (mice_method %in% "norm.boot") {
      m <- 100
    } else {
      m <- 1
    }



    imps <- mice(
      dat_train,
      m = m,
      formulas = formulas,
      method = mice_method,
      printFlag = FALSE,
      maxit = 5
    )



    lmods <- with(imps, lm(y ~ x1 + x2))

    fits <- pool(lmods)
    tidy_fits1 <- tidy(fits, conf.int = TRUE, conf.level = 0.95) %>%
      filter(term == "x1B") %>%
      select(term, estimate, std.error, lower = conf.low, upper = conf.high) %>%
      mutate(true = beta[2], trial = trial) %>%
      mutate(cover = lower < true & true < upper) %>%
      mutate(error = true - estimate) %>%
      mutate(method = method)

    tidy_fits2 <- tidy(fits, conf.int = TRUE, conf.level = 0.95) %>%
      filter(term == "x2") %>%
      select(term, estimate, std.error, lower = conf.low, upper = conf.high) %>%
      mutate(true = beta[3], trial = trial) %>%
      mutate(cover = lower < true & true < upper) %>%
      mutate(error = true - estimate) %>%
      mutate(method = method)

    tidy_fits <- bind_rows(tidy_fits1, tidy_fits2)

    if (mice_method %in% "norm.boot") {
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
    } else {
      preds <- predict(lmods$analyses[[1]], newdata = dat_test, se.fit = TRUE)
      crit <- qt(0.975, df = lmods$analyses[[1]]$df.residual)
      pred_out <- data.frame(qbar = preds$fit)
      se_out <- preds$se.fit
    }


    preds_fits <- data.frame(
      term = "pred",
      estimate = pred_out$qbar,
      std.error = se_out,
      lower = pred_out$qbar - crit * se_out,
      upper = pred_out$qbar + crit * se_out
    )
    preds_fits <- preds_fits %>%
      mutate(true = dat_test$y, trial = trial) %>%
      mutate(cover = lower < true & true < upper) %>%
      mutate(error = true - estimate) %>%
      mutate(method = method)


    bind_rows(tidy_fits, preds_fits)

  }


  pmap(list(
    x = c(TRUE, TRUE, FALSE, FALSE),
    y = c("norm.boot", "norm.predict", "norm.boot", "norm.predict")),
    \(x, y) get_single_trial(trial, x, y)) %>%
    bind_rows()
}



nsim <- 2004
nc <- 50
cl <- makeCluster(nc)
clusterEvalQ(cl, {
  library(tidyverse)
  library(mice)
})
out <- parLapply(cl, 1:nsim, safely(get_trial))
stopCluster(cl)

errors <- map_lgl(out, ~ !is.null(.x$error))
sum(errors)

out_success <- map(out, ~ .x$result)[!errors]

out <- bind_rows(out_success)
out_summary <- out %>%
  group_by(term, method) %>%
  summarize(mbias = mean(error), rmse = sqrt(mean(error^2)), sterr = mean(std.error), cover = mean(cover), .groups = "drop")
out_summary
write_csv(out_summary, here("inst", "output", "simulation", "include_y_summary.csv"))
