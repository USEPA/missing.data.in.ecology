library(tidyverse)
library(here)
library(xtable)

# simulation study

sim_fixed <- read_csv(here("inst", "output", "simulation", "fixed_summary.csv"))
sim_pred <- read_csv(here("inst", "output", "simulation", "pred_summary.csv"))

# data wrangle
sim_pred <- sim_pred %>%
  mutate(term = "pred") %>%
  rename(mse = mspe)

methods_keep <- c("bayesian", "cc", "mean", "norm.nob", "norm.nob_multiple", "norm.boot",
                  "norm.boot_multiple", "norm.predict",
                  "pmm_1D_type0", "pmm_1D_type1", "pmm_1D_type1_multiple",
                  "pmm_3D_type0", "pmm_3D_type0_multiple", "pmm_3D_type1", "pmm_3D_type1_multiple")
sim <- bind_rows(sim_fixed, sim_pred) %>%
  filter(term != "(Intercept)", method %in% methods_keep) %>%
  mutate(rmse = sqrt(mse), SE = sqrt(mvar)) %>%
  select(term, Method = method, Bias = mbias, RMSE = rmse, SE, Cover95 = cover) %>%
  pivot_wider(names_from = term, values_from = -c(term, Method)) %>%
  mutate(Method = case_when(
    Method == "cc" ~ "CCA",
    Method == "mean" ~ "Mean",
    Method == "norm.predict" ~ "Reg",
    Method == "norm.nob" ~ "StReg-S",
    Method == "norm.nob_multiple" ~ "StReg-M",
    Method == "norm.boot" ~ "Boot-S",
    Method == "norm.boot_multiple" ~ "Boot-M",
    Method == "pmm_1D_type0" ~ "PMM-1D-T0",
    Method == "pmm_1D_type1" ~ "PMM-1D-S-T1",
    Method == "pmm_1D_type1_multiple" ~ "PMM-1D-M-T1",
    Method == "pmm_3D_type0" ~ "PMM-3D-S-T0",
    Method == "pmm_3D_type0_multiple" ~ "PMM-3D-M-T0",
    Method == "pmm_3D_type1" ~ "PMM-3D-S-T1",
    Method == "pmm_3D_type1_multiple" ~ "PMM-3D-M-T1",
    Method == "bayesian" ~ "FBDA",
    .default = Method
  )) %>%
  slice(2, 3, 8, 6:7, 4:5, 9:15, 1)


sim_fixed <- sim %>%
  select(Method, Bias_x1B, RMSE_x1B, SE_x1B, Cover95_x1B, Bias_x2, RMSE_x2, SE_x2, Cover95_x2)

sim_pred <- sim %>%
  select(Method, Bias_pred, RMSPE_pred = RMSE_pred, SE_pred, Cover95_pred)

## tables
print(xtable(sim_fixed), include.rownames = FALSE)
print(xtable(sim_pred), include.rownames = FALSE)

## imputation study
sim_imp <- read_csv(here("inst", "output", "simulation", "imp_stat_summary.csv"))

sim_imp <-  sim_imp %>%
  filter(method %in% methods_keep) %>%
  mutate(Method = method) %>%
  select(-method) %>%
  mutate(Method = case_when(
    Method == "cc" ~ "CCA",
    Method == "mean" ~ "Mean",
    Method == "norm.predict" ~ "Reg",
    Method == "norm.nob" ~ "StReg-S",
    Method == "norm.nob_multiple" ~ "StReg-M",
    Method == "norm.boot" ~ "Boot-S",
    Method == "norm.boot_multiple" ~ "Boot-M",
    Method == "pmm_1D_type0" ~ "PMM-1D-T0",
    Method == "pmm_1D_type1" ~ "PMM-1D-S-T1",
    Method == "pmm_1D_type1_multiple" ~ "PMM-1D-M-T1",
    Method == "pmm_3D_type0" ~ "PMM-3D-S-T0",
    Method == "pmm_3D_type0_multiple" ~ "PMM-3D-M-T0",
    Method == "pmm_3D_type1" ~ "PMM-3D-S-T1",
    Method == "pmm_3D_type1_multiple" ~ "PMM-3D-M-T1",
    Method == "bayesian" ~ "FBDA",
    .default = Method
)) %>%
  slice(2, 3, 8, 6:7, 4:5, 9:15, 1) %>%
  mutate(x2_rmse = sqrt(x2_mspe), y_rmspe = sqrt(y_mspe)) %>%
  select(Method, Bias_x2B_imp = x2_bias, RMSE_x2B_imp = x2_rmse, Bias_pred_imp = y_bias, RMSPE_pred_imp = y_rmspe)

sim_imp_fixed <- sim_imp %>%
  left_join(sim_fixed) %>%
  select(Method, contains("x2"))

sim_imp_pred <- sim_imp %>%
  left_join(sim_pred) %>%
  select(Method, contains("pred"))

sim_imp_fixed <- sim_imp_fixed %>%
  relocate(Method, Bias_x2, RMSE_x2, Cover95_x2)
sim_imp_pred <- sim_imp_pred %>%
  relocate(Method, Bias_pred, RMSPE_pred, Cover95_pred)
print(xtable(sim_imp_fixed), include.rownames = FALSE)
print(xtable(sim_imp_pred), include.rownames = FALSE)

# water present response

water1 <- read_csv(here("inst", "output", "water-present", "water_present.csv")) %>%
  mutate(Parameter = c("Water Present", "Water Not Present"))
water2 <- read_csv(here("inst", "output", "water-present", "ntl_cond_water_present.csv")) %>%
  mutate(Parameter = str_c("Nitrogen ",Category, ", Water Present"))
water <- bind_rows(water1, water2) %>%
  mutate(across(c(Estimate.P, StdError.P, LCB95Pct.P, UCB95Pct.P), \(x) format(round(x, digits = 2)))) %>%
  mutate(across(c(Estimate.P, StdError.P, LCB95Pct.P, UCB95Pct.P), \(x) str_c(x, "%"))) %>%
  mutate("Estimate (SE)" = str_c(Estimate.P, " (", StdError.P, ")")) %>%
  select(Parameter, `Estimate (SE)`, LCB95 = LCB95Pct.P, UCB95 = UCB95Pct.P)

## tables
print(xtable(water), include.rownames = FALSE)


# water present explanatory
cvmmi <- read_csv(here("inst", "output", "vmmi", "fixed_contingency.csv")) %>%
  filter(term != "(Intercept)") %>%
  mutate(across(c(estimate:conf.high), \(x) format(round(x, digits = 2)))) %>%
  mutate("Estimate (SE)" = str_c(estimate, " (", std.error, ")")) %>%
  select(Parameter = term, `Estimate (SE)`, LCB95 = conf.low, UCB95 = conf.high, `p-value` = p.value)

print(xtable(cvmmi), include.rownames = FALSE)

#
fixed <- read_csv(here("inst", "output", "vmmi", "fixed.csv")) %>%
  mutate(conf.low = estimate - qnorm(0.975) * std.error, conf.high = estimate + qnorm(0.975) * std.error) %>%
  filter(term != "(Intercept)") %>%
  relocate(type) %>%
  mutate(across(c(estimate:conf.high), \(x) format(round(x, digits = 2)))) %>%
  mutate("Estimate (SE)" = str_c(estimate, " (", std.error, ")")) %>%
  select(type, term, `Estimate (SE)`) %>%
  pivot_wider(names_from = type, values_from = c(`Estimate (SE)`)) %>%
  select(Parameter = term, `none-all`, `none-missing`, `none-cc`, `sp-all`, `sp-cc`, `sp-missing`)

print(xtable(fixed %>% select(Parameter, `none-all`, `none-missing`, `none-cc`)), include.rownames = FALSE)
print(xtable(fixed %>% select(Parameter, `sp-all`, `sp-missing`, `sp-cc`)), include.rownames = FALSE)

stats <- read_csv(here("inst", "output", "vmmi", "stats.csv")) %>%
  mutate(Model = if_else(str_detect(type, "sp"), "Spatial", "Nonspatial")) %>%
  mutate(Approach = case_when(
    str_detect(type, "missing") ~ "Multiple Imputation",
    str_detect(type, "all") ~ "All Data",
    str_detect(type, "cc") ~ "Complete Case"
  )) %>%
  mutate(across(c(bias:BIC), \(x) format(round(x, digits = 2)))) %>%
  select(Model, Approach, Bias = bias, RMSPE, R2 = cor2) %>%
  mutate(Approach = factor(Approach, levels = c("All Data", "Multiple Imputation", "Complete Case"))) %>%
  arrange(Model, Approach)

print(xtable(stats), include.rownames = FALSE)


# vmmi2
fixed <- read_csv(here("inst", "output", "vmmi", "fixed2.csv")) %>%
  mutate(conf.low = estimate - qnorm(0.975) * std.error, conf.high = estimate + qnorm(0.975) * std.error) %>%
  filter(term != "(Intercept)") %>%
  relocate(type) %>%
  mutate(across(c(estimate:conf.high), \(x) format(round(x, digits = 3)))) %>%
  mutate("Estimate (SE)" = str_c(estimate, " (", std.error, ")")) %>%
  select(type, term, `Estimate (SE)`) %>%
  pivot_wider(names_from = type, values_from = c(`Estimate (SE)`)) %>%
  select(Parameter = term, `none-all`, `none-missing`, `none-cc`, `sp-all`, `sp-cc`, `sp-missing`)


print(xtable(fixed %>% select(Parameter, `sp-all`, `sp-missing`, `sp-cc`)), include.rownames = FALSE)

stats <- read_csv(here("inst", "output", "vmmi", "stats2.csv"))
print(xtable(stats %>% select(type, MBias = bias, RMSPE, PCover95 = cover95), digits = 3), include.rownames = FALSE)


# crossing output
crossing <- read_csv(here("inst", "output", "simulation", "crossing_summary.csv")) %>%
  arrange(desc(term)) %>%
  separate(method, c("x2_imp", "x2_mod", "y_imp"), sep = "_") %>%
  mutate(across(c(x2_imp, x2_mod, y_imp), \(x) map_chr(str_split(x, "-"), \(y) y[3]))) %>%
  mutate(across(c(x2_imp, x2_mod, y_imp), str_to_title)) %>%
  pivot_wider(names_from = term, values_from = c(mbias, rmse, sterr, cover), names_vary = "slowest") %>%
  relocate(x2_imp, y_imp, x2_mod) %>%
  arrange(desc(x2_imp), desc(y_imp), desc(x2_mod))

print(xtable(crossing), include.rownames = FALSE, digits = 3)


# y output (fixed effects)
include_y <- read_csv(here("inst", "output", "simulation", "include_y_summary.csv")) %>%
  filter(term %in% c("x1B", "x2")) %>%
  separate(method, c("method", "include-y"), sep = "_") %>%
  mutate(method = if_else(method == "norm.boot", "Boot-M", "Reg")) %>%
  pivot_wider(names_from = term, values_from = c(mbias, rmse, sterr, cover), names_vary = "slowest") %>%
  arrange(desc(method))

print(xtable(include_y), include.rownames = FALSE, digits = 3)

# y output (prediction)
include_y <- read_csv(here("inst", "output", "simulation", "include_y_summary.csv")) %>%
  filter(term %in% c("pred")) %>%
  arrange(desc(term)) %>%
  separate(method, c("method", "include-y"), sep = "_") %>%
  mutate(method = if_else(method == "norm.boot", "Boot-M", "Reg")) %>%
  select(-term) %>%
  arrange(desc(method))

print(xtable(include_y), include.rownames = FALSE, digits = 3)


# simulation for supporting information

sim_fixed <- read_csv(here("inst", "output", "simulation", "fixed_summary.csv"))
sim_pred <- read_csv(here("inst", "output", "simulation", "pred_summary.csv"))

# data wrangle
sim_pred <- sim_pred %>%
  mutate(term = "pred") %>%
  rename(mse = mspe)

methods_keep <- c("bayesian", "cc", "mean", "norm.nob", "norm.nob_multiple", "norm.boot",
                  "norm.boot_multiple", "norm.predict",
                  "pmm_1D_type0", "rf", "rf_multiple", "lasso.norm", "lasso.norm_multiple",
                  "norm", "norm_multiple",
                  "pmm_3D_type0", "pmm_3D_type0_multiple", "pmm_3D_type1", "pmm_3D_type1_multiple")
sim <- bind_rows(sim_fixed, sim_pred) %>%
  filter(term != "(Intercept)", method %in% methods_keep) %>%
  mutate(rmse = sqrt(mse)) %>%
  select(term, Method = method, Bias = mbias, RMSE = rmse, Cover95 = cover) %>%
  pivot_wider(names_from = term, values_from = -c(term, Method)) %>%
  mutate(Method = case_when(
    Method == "cc" ~ "CCA",
    Method == "mean" ~ "Mean",
    Method == "norm.predict" ~ "Reg",
    Method == "norm.nob" ~ "StReg-S",
    Method == "norm.nob_multiple" ~ "StReg-M",
    Method == "norm.boot" ~ "Boot-S",
    Method == "norm.boot_multiple" ~ "Boot-M",
    Method == "norm" ~ "BReg-S",
    Method == "lasso.norm" ~ "LReg-S",
    Method == "rf" ~ "RF-S",
    Method == "pmm_1D_type0" ~ "PMM-1D-T0",
    Method == "pmm_1D_type1" ~ "PMM-1D-S-T1",
    Method == "pmm_1D_type1_multiple" ~ "PMM-1D-M-T1",
    Method == "pmm_3D_type0" ~ "PMM-3D-S-T0",
    Method == "pmm_3D_type0_multiple" ~ "PMM-3D-M-T0",
    Method == "pmm_3D_type1" ~ "PMM-3D-S-T1",
    Method == "pmm_3D_type1_multiple" ~ "PMM-3D-M-T1",
    Method == "norm_multiple" ~ "BReg-M",
    Method == "lasso.norm_multiple" ~ "LReg-M",
    Method == "rf_multiple" ~ "RF-M",
    Method == "bayesian" ~ "FBDA",
    .default = Method
  )) %>%
  slice(2, 5, 11, 13, 9, 14, 7, 16, 6, 3, 18, 10, 15, 8, 17, 12, 4, 19, 1)


sim_fixed <- sim %>%
  select(Method, Bias_x1B, RMSE_x1B, Cover95_x1B, Bias_x2, RMSE_x2, Cover95_x2)

sim_pred <- sim %>%
  select(Method, Bias_pred, RMSPE_pred = RMSE_pred, Cover95_pred)

## tables
print(xtable(sim_fixed), include.rownames = FALSE)
print(xtable(sim_pred), include.rownames = FALSE)

# MNAR output
mnar_small <- read_csv(here("inst", "output", "simulation", "mnar_summary_small.csv")) %>%
  mutate(type = "small")
mnar_moderate <- read_csv(here("inst", "output", "simulation", "mnar_summary_moderate.csv")) %>%
  mutate(type = "moderate")
mnar_large <- read_csv(here("inst", "output", "simulation", "mnar_summary_large.csv")) %>%
  mutate(type = "large")
mnar <- bind_rows(mnar_small, mnar_moderate, mnar_large) %>%
  arrange(type) %>%
  mutate(model = case_when(
    model == "MAR/MCAR" ~ "MAR",
    model == "heckman" ~ "Heckman",
    .default = model
  )) %>%
  mutate(type = str_to_title(type))
mnar_fixed <- mnar %>%
  filter(variable != "y_ho", variable != "marg_b_y[1]")
mnar_pred <- mnar %>%
  filter(variable == "y_ho")

mnar_fixed %>%
  filter(model %in% c("MAR", "Heckman")) %>%
  pivot_wider(names_from = variable, values_from = c(mbias, rmse, se, cover95), names_vary = "slowest") %>%
  xtable() %>%
  print(include.rownames = FALSE)

mnar_pred %>%
  filter(model %in% c("MAR", "Heckman")) %>%
  pivot_wider(names_from = variable, values_from = c(mbias, rmse, se, cover95), names_vary = "slowest") %>%
  xtable() %>%
  print(include.rownames = FALSE)

mnar_fixed %>%
  pivot_wider(names_from = variable, values_from = c(mbias, rmse, se, cover95), names_vary = "slowest") %>%
  xtable() %>%
  print(include.rownames = FALSE)

mnar_pred %>%
  pivot_wider(names_from = variable, values_from = c(mbias, rmse, se, cover95), names_vary = "slowest") %>%
  xtable() %>%
  print(include.rownames = FALSE)
