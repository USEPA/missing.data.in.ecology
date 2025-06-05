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
  select(Method, Bias_x1B, RMSE_x1B, Cover95_x1B, Bias_x2, RMSE_x2, Cover95_x2)

sim_pred <- sim %>%
  select(Method, Bias_pred, RMSPE_pred = RMSE_pred, Cover95_pred)

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
  pivot_wider(names_from = term, values_from = c(mbias, rmse, sterr, cover), names_vary = "slowest")

print(xtable(crossing), include.rownames = FALSE)
