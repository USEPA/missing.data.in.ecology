library(tidyverse)
library(here)
library(mice)
library(sf)
library(spmodel)
library(emmeans)
library(parallel)

# set initial seed
set.seed(10)

#### initial functions
##### impose missing data
sim_missing <- function(IS_WOODY) {

  if (IS_WOODY) {
    p <- 0.3
  }

  if (!IS_WOODY) {
    p <- 0.6
  }

  binvar <- rbinom(1, 1, p)

}

get_imputed <- function(x, formula, imps, spcov_type, dat_train, dat_test) {
  dat_train <- complete(imps, x)
  mod <- splm(formula, data = dat_train, spcov_type = spcov_type, xcoord = XCOORD, ycoord = YCOORD)
  pred <- predict(mod, newdata = dat_test, se.fit = TRUE)
  list(mod = mod, coefs = coef(mod), vcoefs = diag(vcov(mod)), pred = pred$fit, pred_se = pred$se.fit, true_pred = dat_test[["VMMI_2016"]])
}

get_params <- function(index, fits) {
  name <- names(fits[[1]]$coefs)[index]
  coefs <- sapply(1:m, function(x) fits[[x]]$coefs[index])
  vcoefs <- sapply(1:m, function(x) fits[[x]]$vcoefs[index])
  pools <- pool.scalar(coefs, vcoefs)
  dat <- data.frame(term = name, estimate = pools$qbar, std.error = sqrt(pools$t))
  dat$statistic <- dat$estimate / dat$std.error
  dat$p.value <- 2 * pt(abs(dat$statistic), df = Inf, lower.tail = FALSE)
  dat
}

get_stats <- function(fits, p) {
  ests <- tibble::as_tibble(map_dfr(seq_len(p), get_params, fits = fits))
  pred_fits <- bind_cols(map(fits, \(x) x$pred))
  n_pred <- NROW(pred_fits)
  pred_vars <- bind_cols(map(fits, \(x) x$pred_se^2))
  pred_out <- map(seq_len(n_pred), \(x) pool.scalar(Q = unlist(pred_fits[x, ]), U = unlist(pred_vars[x, ]), n = n_pred))
  pred_fits <- tibble::tibble(
    VMMI_2016 = fits[[1]]$true_pred,
    fit = map_dbl(pred_out, \(x) x$qbar),
    resid = VMMI_2016 - fit,
    se = sqrt(map_dbl(pred_out, \(x) x$t)),
    lower95 = fit - qnorm(0.975) * se,
    upper95 = fit + qnorm(0.975) * se,
    cover = if_else(between(VMMI_2016, lower95, upper95), 1, 0)
  )
  fit_stat <- tibble::tibble(
    bias = mean(pred_fits$resid),
    MSPE = mean(pred_fits$resid^2),
    RMSPE = sqrt(MSPE),
    cor = cor(pred_fits$VMMI_2016, pred_fits$fit),
    cover95 = mean(pred_fits$cover)
  )
  list(fixed = ests, fit_stat = fit_stat, pred_fits = pred_fits)
}

######################### start


dat_full <- read_csv(here("inst", "data", "nwca_2016.csv"), guess_max = Inf)

dat <- dat_full %>%
  mutate(SAMPLEABLE_WATER = if_else(SAMPLEABLE_H2O == "Yes", 1, 0)) %>%
  mutate(IS_WOODY = if_else(WETCLS_GRP %in% c("PRLW", "EW"), 1, 0)) %>%
  mutate(NTL_RESULT = if_else(is.na(NTL_RESULT), 0, log(NTL_RESULT))) %>%
  select(VMMI_2016, IS_WOODY, SAMPLEABLE_WATER, NTL_RESULT,
         STRESS_VEGRMV, PALT_SOHARD, PALT_SOMODF, LON_DD83, LAT_DD83, XCOORD, YCOORD) %>%
  mutate(PALT_SOMODF = as.numeric(PALT_SOMODF), PALT_SOHARD = as.numeric(PALT_SOHARD),
         STRESS_VEGRMV = factor(STRESS_VEGRMV, levels = c("Low", "Moderate", "High"))) %>%
  filter(!is.na(PALT_SOMODF), !is.na(VMMI_2016))

dat$PALT_SOMODF_PRESENT <- if_else(dat$PALT_SOMODF > 0, 1, 0)

train_prop <- 2/3
n_dat <- NROW(dat)
train_index <- sample(seq(1, n_dat), round(train_prop * n_dat))
dat_train <- dat[train_index, ]
dat_test <- dat[-train_index, ]


# fuller model
form <- VMMI_2016 ~ PALT_SOHARD + IS_WOODY

# linear model
lmod <- splm(form, data = dat_train, spcov_type = "none")
summary(lmod)
loocv(lmod)

lmod_preds <- predict(lmod, newdata = dat_test, interval = "prediction")


# spatial model
spmod <- splm(form, data = dat_train, spcov_type = "exponential", xcoord = XCOORD, ycoord = YCOORD)
summary(spmod)
loocv(spmod)
spmod_preds <- predict(spmod, newdata = dat_test, interval = "prediction")


# compare models
glances(lmod, spmod)



dat_train <- dat_train %>%
  mutate(
    PALT_SOHARD_MISS = map_dbl(IS_WOODY, sim_missing)# ,
    # VMMI_2016_MISS = map_dbl(IS_WOODY, sim_missing)
  ) %>%
  mutate(
    PALT_SOHARD = if_else(PALT_SOHARD_MISS == 1, NA, PALT_SOHARD)# ,
    # VMMI_2016 = if_else(VMMI_2016_MISS == 1, NA, VMMI_2016)
  )

lmod_cc <- splm(form, data = drop_na(dat_train), spcov_type = "none")
summary(lmod_cc)
loocv(lmod_cc)
lmod_cc_preds <- predict(lmod_cc, newdata = dat_test, interval = "prediction")

spmod_cc <- splm(form, data = drop_na(dat_train), spcov_type = "exponential", xcoord = XCOORD, ycoord = YCOORD)
summary(spmod_cc)
loocv(spmod_cc)
spmod_cc_preds <- predict(spmod_cc, newdata = dat_test, interval = "prediction")


m <- 100
formulas <- list(
  # PALT_SOHARD ~ VMMI_2016
  PALT_SOHARD ~ IS_WOODY + VMMI_2016# ,
  # VMMI_2016 ~ PALT_SOHARD + IS_WOODY
)
formulas <- name.formulas(formulas)

imps <- mice(
  dat_train,
  m = m,
  formulas = formulas,
  method = "norm.boot", # or pmm
  printFlag = FALSE
)

nc <- 50
cl <- makeCluster(nc)
clusterEvalQ(cl, {
  library(spmodel)
  library(mice)
})
fits_sp <- parLapply(cl, 1:m, get_imputed, formula = form, imps = imps, spcov_type = "exponential", dat_train = dat_train, dat_test = dat_test)
fits_none <- parLapply(cl, 1:m, get_imputed, formula = form, imps = imps, spcov_type = "none", dat_train = dat_train, dat_test = dat_test)
stopCluster(cl)

p <- length(tidy(spmod)$term)
stats_sp <- get_stats(fits_sp, p)
stats_none <- get_stats(fits_none, p)

fixed_out <- bind_rows(
  stats_sp$fixed %>% mutate(type = "sp-missing"),
  tidy(spmod) %>% mutate(type = "sp-all"),
  tidy(spmod_cc) %>% mutate(type = "sp-cc"),
  stats_none$fixed %>% mutate(type = "none-missing"),
  tidy(lmod) %>% mutate(type = "none-all"),
  tidy(lmod_cc) %>% mutate(type = "none-cc")
) %>%
  arrange(term)


lmod_pred_df <- tibble(
  bias = mean(lmod_preds[, "fit"] - dat_test$VMMI_2016),
  MSPE = mean((lmod_preds[, "fit"] - dat_test$VMMI_2016)^2),
  RMSPE = sqrt(MSPE),
  cor = cor(lmod_preds[, "fit"], dat_test$VMMI_2016),
  cover95 = mean(if_else(dat_test$VMMI_2016 >= lmod_preds[, "lwr"] & dat_test$VMMI_2016 <= lmod_preds[, "upr"], 1, 0))
)

spmod_pred_df <- tibble(
  bias = mean(spmod_preds[, "fit"] - dat_test$VMMI_2016),
  MSPE = mean((spmod_preds[, "fit"] - dat_test$VMMI_2016)^2),
  RMSPE = sqrt(MSPE),
  cor = cor(spmod_preds[, "fit"], dat_test$VMMI_2016),
  cover95 = mean(if_else(dat_test$VMMI_2016 >= spmod_preds[, "lwr"] & dat_test$VMMI_2016 <= spmod_preds[, "upr"], 1, 0))
)

lmod_cc_pred_df <- tibble(
  bias = mean(lmod_cc_preds[, "fit"] - dat_test$VMMI_2016),
  MSPE = mean((lmod_cc_preds[, "fit"] - dat_test$VMMI_2016)^2),
  RMSPE = sqrt(MSPE),
  cor = cor(lmod_cc_preds[, "fit"], dat_test$VMMI_2016),
  cover95 = mean(if_else(dat_test$VMMI_2016 >= lmod_cc_preds[, "lwr"] & dat_test$VMMI_2016 <= lmod_cc_preds[, "upr"], 1, 0))
)

spmod_cc_pred_df <- tibble(
  bias = mean(spmod_cc_preds[, "fit"] - dat_test$VMMI_2016),
  MSPE = mean((spmod_cc_preds[, "fit"] - dat_test$VMMI_2016)^2),
  RMSPE = sqrt(MSPE),
  cor = cor(spmod_cc_preds[, "fit"], dat_test$VMMI_2016),
  cover95 = mean(if_else(dat_test$VMMI_2016 >= spmod_cc_preds[, "lwr"] & dat_test$VMMI_2016 <= spmod_cc_preds[, "upr"], 1, 0))
)


stats_out <- bind_rows(
  lmod_pred_df %>% mutate(type = "none-all"),
  spmod_pred_df %>% mutate(type = "spmod-all"),
  lmod_cc_pred_df %>% mutate(type = "none-cc"),
  spmod_cc_pred_df %>% mutate(type = "spmod-cc"),
  stats_none$fit_stat %>% mutate(type = "none-missing"),
  stats_sp$fit_stat %>% mutate(type = "sp-missing")
)

stats_out <- stats_out %>%
  mutate(rank_MSPE = rank(MSPE), rank_cor = rank(desc(cor))) %>%
  arrange(rank_MSPE)

lik_out <- tibble(
  AIC = c(AIC(lmod), AIC(spmod), AIC(lmod_cc), AIC(spmod_cc)),
  BIC = c(BIC(lmod), BIC(spmod), BIC(lmod_cc), BIC(spmod_cc)),
  type = c("none-all", "sp-all", "none-cc", "sp-cc")
)

write_csv(fixed_out, str_c(here("inst", "output", "vmmi"), "/fixed2.csv"))
write_csv(stats_out, str_c(here("inst", "output", "vmmi"), "/stats2.csv"))
write_csv(lik_out, str_c(here("inst", "output", "vmmi"), "/lik2.csv"))
lmods <- with(imps, lm(VMMI_2016 ~ IS_WOODY + PALT_SOHARD))
fixed_out_lm <- tidy(pool(lmods))[, 1:5]
write_csv(fixed_out_lm, str_c(here("inst", "output", "vmmi"), "/fixed_lm2.csv"))
