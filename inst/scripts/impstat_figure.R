library(tidyverse)
library(here)

remove_vals <- c("cc", "lasso.norm", "lasso.norm_multiple", "rf", "rf_multiple", "norm", "norm_multiple")

impdat <- read_csv(here("inst", "output", "simulation", "imp_stat_summary.csv")) %>%
  filter(! method %in% remove_vals) %>%
  mutate(imp_RMSPE = sqrt(x2_mspe)) %>%
  select(method, imp_RMSPE)

beta2dat <- read_csv(here("inst", "output", "simulation", "fixed_summary.csv")) %>%
  filter(term == "x2", ! method %in% remove_vals) %>%
  mutate(beta2_RMSE = sqrt(mse)) %>%
  select(method, beta2_RMSE)

dat <- left_join(impdat, beta2dat) %>%
  mutate(Method = case_when(
    method == "mean" ~ "Mean",
    method == "norm.predict" ~ "Reg",
    method == "pmm_1D_type0" ~ "NN",
    .default = "Other"
  )) %>%
  mutate(Method = factor(Method, levels = c("Mean", "Reg", "NN", "Other")))


okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
p1 <- ggplot(dat, aes(x = beta2_RMSE, y = imp_RMSPE, color = Method, shape = Method)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0.35) +
  labs(x = expression(~beta[2]~Parameter~RMSE), y = expression(x[2]~Imputation~RMSE)) +
  scale_color_manual(values = okabe[1:4]) +
  scale_shape_manual(values = c(19, 19, 19, 17)) +
  theme_bw(base_size = 14)

ggsave(
  filename = here("inst", "figures", "imp-vs-model-comp.jpeg"),
  plot = p1,
  dpi = 300,
  height = 4.98,
  width = 6.39
)
