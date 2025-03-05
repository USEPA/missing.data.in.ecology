library(tidyverse)
library(here)
library(mice)
library(sf)
library(spmodel)
library(emmeans)
library(parallel)

# set initial seed
set.seed(0)

#### initial functions
##### impose missing data
sim_missing <- function(IS_WOODY, SAMPLEABLE_WATER) {

  if (IS_WOODY & SAMPLEABLE_WATER) {
    p <- 0.2
  }

  if (!IS_WOODY & SAMPLEABLE_WATER) {
    p <- 0.4
  }

  if (IS_WOODY & !SAMPLEABLE_WATER) {
    p <- 0.6
  }

  if (!IS_WOODY & !SAMPLEABLE_WATER) {
    p <- 0.8
  }

  binvar <- rbinom(1, 1, p)

}

get_imputed <- function(x, formula, imps, spcov_type) {
  dat <- complete(imps, x)
  mod <- splm(formula, data = dat, spcov_type = spcov_type, xcoord = XCOORD, ycoord = YCOORD)
  loocv_fit <- loocv(mod, cv_predict = TRUE, se.fit = TRUE)
  cv <- loocv_fit$cv_predict
  cv_var <- loocv_fit$se.fit^2
  list(mod = mod, coefs = coef(mod), vcoefs = diag(vcov(mod)), cv = cv, cv_var = cv_var, AIC = AIC(mod), BIC = BIC(mod))
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

get_stats <- function(fits, data, p) {
  ests <- tibble::as_tibble(map_dfr(seq_len(p), get_params, fits = fits))
  n <- fits[[1]]$mod$n
  loocv_fits <- bind_cols(map(fits, \(x) x$cv))
  loocv_vars <- bind_cols(map(fits, \(x) x$cv_var))
  loocv_out <- map(seq_len(n), \(x) pool.scalar(Q = unlist(loocv_fits[x, ]), U = unlist(loocv_vars[x, ]), n = n))
  loocv_fits <- tibble::tibble(
    VMMI_2016 = data$VMMI_2016,
    fit = map_dbl(loocv_out, \(x) x$qbar),
    resid = VMMI_2016 - fit,
    se = sqrt(map_dbl(loocv_out, \(x) x$t)),
    lower95 = fit - qnorm(0.975) * se,
    upper95 = fit + qnorm(0.975) * se,
    cover = if_else(between(VMMI_2016, lower95, upper95), 1, 0)
  )
  fit_stat <- tibble::tibble(
    bias = mean(loocv_fits$resid),
    MSPE = mean(loocv_fits$resid^2),
    RMSPE = sqrt(MSPE),
    cor2 = cor(loocv_fits$VMMI_2016, loocv_fits$fit)^2,
    cover95 = mean(loocv_fits$cover),
    AIC = mean(map_dbl(fits, \(x) x$AIC)),
    BIC = mean(map_dbl(fits, \(x) x$BIC))
  )
  list(fixed = ests, fit_stat = fit_stat, loocv_fits = loocv_fits)
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

# linear model fit to all data; contingency filter
# lmod <- splm(VMMI_2016 ~ SAMPLEABLE_WATER + SAMPLEABLE_WATER:NTL_RESULT,
#              data = dat, spcov_type = "none")
# summary(lmod)

dat$PALT_SOMODF_PRESENT <- if_else(dat$PALT_SOMODF > 0, 1, 0)
lmod <- splm(VMMI_2016 ~ PALT_SOMODF_PRESENT + SAMPLEABLE_WATER + SAMPLEABLE_WATER:NTL_RESULT ,
             data = dat, spcov_type = "none")
summary(lmod)
write_csv(tidy(lmod, conf.int = TRUE), str_c(here("inst", "output", "vmmi"), "/fixed_contingency.csv"))
# new
# lmod <- splm(VMMI_2016 ~ SAMPLEABLE_WATER:NTL_RESULT + PALT_SOMODF_PRESENT,
#              data = dat, spcov_type = "none")
# summary(lmod)
# lmod <- splm(VMMI_2016 ~ SAMPLEABLE_WATER:NTL_RESULT + SAMPLEABLE_WATER + PALT_SOMODF_PRESENT,
#              data = dat, spcov_type = "none")
# summary(lmod)

# fuller model
form <- VMMI_2016 ~ IS_WOODY + SAMPLEABLE_WATER + SAMPLEABLE_WATER:NTL_RESULT +
  PALT_SOHARD + PALT_SOMODF_PRESENT + STRESS_VEGRMV

# linear model
lmod <- splm(form, data = dat, spcov_type = "none")
summary(lmod)
loocv(lmod)

# spatial model
spmod <- splm(form, data = dat, spcov_type = "exponential", xcoord = XCOORD, ycoord = YCOORD)
summary(spmod)
loocv(spmod)

# compare models
glances(lmod, spmod)


# add missing indicator
dat <- dat %>%
  mutate(
    SOHARD_MISS = map2_dbl(IS_WOODY, SAMPLEABLE_WATER, sim_missing),
    SOMODF_MISS = map2_dbl(IS_WOODY, SAMPLEABLE_WATER, sim_missing),
    VEGRMV_MISS = map2_dbl(IS_WOODY, SAMPLEABLE_WATER, sim_missing),
  ) %>%
  mutate(
    PALT_SOHARD = if_else(SOHARD_MISS == 1, NA, PALT_SOHARD),
    PALT_SOMODF_PRESENT = if_else(SOMODF_MISS == 1, NA, PALT_SOMODF_PRESENT),
    STRESS_VEGRMV = factor(if_else(VEGRMV_MISS == 1, NA, STRESS_VEGRMV), levels = c("Low", "Moderate", "High", NA))
  )

dat %>%
  group_by(IS_WOODY, SAMPLEABLE_WATER) %>%
  summarize(
    PALT_SOHARD = mean(is.na(PALT_SOHARD)),
    PALT_SOMODF_PRESENT = mean(is.na(PALT_SOMODF_PRESENT)),
    STRESS_VEGRMV = mean(is.na(STRESS_VEGRMV))
  )
dat <- dat %>%
  select(VMMI_2016, IS_WOODY, SAMPLEABLE_WATER, NTL_RESULT, PALT_SOMODF_PRESENT, PALT_SOHARD,
         STRESS_VEGRMV, XCOORD, YCOORD)

# linear model
lmod2 <- splm(form, data = drop_na(dat), spcov_type = "none")
summary(lmod2)
loocv(lmod2)

# spatial mod2el
spmod2 <- splm(form, data = drop_na(dat), spcov_type = "exponential", xcoord = XCOORD, ycoord = YCOORD)
summary(spmod2)
loocv(spmod2)

# compare mod2els
glances(lmod2, spmod2)

# categorical response need to be a factor
m <- 100
formulas <- list(
  PALT_SOHARD ~ IS_WOODY + SAMPLEABLE_WATER + SAMPLEABLE_WATER:NTL_RESULT + PALT_SOMODF_PRESENT + STRESS_VEGRMV + VMMI_2016,
  PALT_SOMODF_PRESENT ~ IS_WOODY + SAMPLEABLE_WATER + SAMPLEABLE_WATER:NTL_RESULT + PALT_SOHARD + STRESS_VEGRMV + VMMI_2016,
  STRESS_VEGRMV ~ IS_WOODY + SAMPLEABLE_WATER + SAMPLEABLE_WATER:NTL_RESULT + PALT_SOHARD + PALT_SOMODF_PRESENT + VMMI_2016
)
formulas <- name.formulas(formulas)

imps <- mice(
  dat,
  m = m,
  formulas = formulas,
  method = "pmm",
  printFlag = FALSE
)

nc <- 10
cl <- makeCluster(nc)
clusterEvalQ(cl, {
  library(spmodel)
  library(mice)
})
fits_sp <- parLapply(cl, 1:m, get_imputed, formula = form, imps = imps, spcov_type = "exponential")
fits_none <- parLapply(cl, 1:m, get_imputed, formula = form, imps = imps, spcov_type = "none")
stopCluster(cl)

p <- length(tidy(spmod)$term)
stats_sp <- get_stats(fits_sp, dat, p)
stats_none <- get_stats(fits_none, dat, p)

fixed_out <- bind_rows(
  stats_sp$fixed %>% mutate(type = "sp-missing"),
  tidy(spmod) %>% mutate(type = "sp-all"),
  tidy(spmod2) %>% mutate(type = "sp-cc"),
  stats_none$fixed %>% mutate(type = "none-missing"),
  tidy(lmod) %>% mutate(type = "none-all"),
  tidy(lmod2) %>% mutate(type = "none-cc")
) %>%
  arrange(term)

loocv_spmod_all <- loocv(spmod, cv_predict = TRUE, se.fit = TRUE)
loocv_spmod_all$stats$cover95 <- mean(between(
  spmod$obdata$VMMI_2016,
  loocv_spmod_all$cv_predict - qnorm(0.975) * loocv_spmod_all$se.fit,
  loocv_spmod_all$cv_predict + qnorm(0.975) * loocv_spmod_all$se.fit
))
loocv_spmod2_all <- loocv(spmod2, cv_predict = TRUE, se.fit = TRUE)
loocv_spmod2_all$stats$cover95 <- mean(between(
  spmod2$obdata$VMMI_2016,
  loocv_spmod2_all$cv_predict - qnorm(0.975) * loocv_spmod2_all$se.fit,
  loocv_spmod2_all$cv_predict + qnorm(0.975) * loocv_spmod2_all$se.fit
))
loocv_lmod_all <- loocv(lmod, cv_predict = TRUE, se.fit = TRUE)
loocv_lmod_all$stats$cover95 <- mean(between(
  lmod$obdata$VMMI_2016,
  loocv_lmod_all$cv_predict - qnorm(0.975) * loocv_lmod_all$se.fit,
  loocv_lmod_all$cv_predict + qnorm(0.975) * loocv_lmod_all$se.fit
))
loocv_lmod2_all <- loocv(lmod2, cv_predict = TRUE, se.fit = TRUE)
loocv_lmod2_all$stats$cover95 <- mean(between(
  lmod2$obdata$VMMI_2016,
  loocv_lmod2_all$cv_predict - qnorm(0.975) * loocv_lmod2_all$se.fit,
  loocv_lmod2_all$cv_predict + qnorm(0.975) * loocv_lmod2_all$se.fit
))
stats_out <- bind_rows(
  stats_sp$fit_stat %>% mutate(type = "sp-missing"),
  loocv_spmod_all$stats %>% mutate(AIC = AIC(spmod), BIC = BIC(spmod), type = "sp-all"),
  loocv_spmod2_all$stats %>% mutate(AIC = AIC(spmod2), BIC = BIC(spmod2), type = "sp-cc"),
  stats_none$fit_stat %>% mutate(type = "none-missing"),
  loocv_lmod_all$stats %>% mutate(AIC = AIC(lmod), BIC = BIC(lmod), type = "none-all"),
  loocv_lmod2_all$stats %>% mutate(AIC = AIC(lmod2), BIC = BIC(lmod2), type = "none-cc")
)

loocv_out_all <- tibble(
  response = dat$VMMI_2016,
  "sp-all" = loocv(spmod, cv_predict = TRUE)$cv_predict,
  "none-all" = loocv(lmod, cv_predict = TRUE)$cv_predict,
)

loocv_out_miss <- tibble(
  response = dat$VMMI_2016,
  "sp-missing" = stats_sp$loocv_fits$fit,
  "none-missing" = stats_none$loocv_fits$fit,
)

loocv_out_cc <- tibble(
  response = drop_na(dat)$VMMI_2016,
  "sp-cc" = loocv(spmod2, cv_predict = TRUE)$cv_predict,
  "none-cc" = loocv(lmod2, cv_predict = TRUE)$cv_predict
)

write_csv(tidy(lmod, conf.int = TRUE), str_c(here("inst", "output", "vmmi"), "/fixed_contingency.csv"))
write_csv(fixed_out, str_c(here("inst", "output", "vmmi"), "/fixed.csv"))
write_csv(stats_out, str_c(here("inst", "output", "vmmi"), "/stats.csv"))
# loocv
# write_csv(loocv_out_all, str_c(here("inst", "output", "vmmi"), "/loocv_all.csv"))
# write_csv(loocv_out_miss, str_c(here("inst", "output", "vmmi"), "/loocv_miss.csv"))
# write_csv(loocv_out_cc, str_c(here("inst", "output", "vmmi"), "/loocv_cc.csv"))

# check none vs this
lmods <- with(imps, lm(VMMI_2016 ~ IS_WOODY + SAMPLEABLE_WATER + SAMPLEABLE_WATER:NTL_RESULT +
                         PALT_SOHARD + PALT_SOMODF_PRESENT + STRESS_VEGRMV))
fixed_out_lm <- tidy(pool(lmods))[, 1:5]
write_csv(fixed_out_lm, str_c(here("inst", "output", "vmmi"), "/fixed_lm.csv"))

fixed_out
stats_out
