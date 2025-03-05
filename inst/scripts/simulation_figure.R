library(tidyverse)
library(mice)
library(brms)

# make simulation figure


set.seed(17)
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
  mutate(x2_BINOM = map_dbl(x1, get_binom, p1 = 0.35, p2 = 0.70)) %>%
  mutate(y_BINOM = map_dbl(x1, get_binom, p1 = 0.70, p2 = 0.35)) %>%
  mutate(x2 = if_else(x2_BINOM == 1, NA, x2)) %>%
  mutate(y = if_else(y_BINOM == 1, NA, y))

formulas <- list(
  y ~ x1 + x2,
  x2 ~ x1 + y
)

formulas <- name.formulas(formulas)

dat_plot <- dat_train %>%
  mutate(Status = factor(if_else(is.na(x2), "Imputed", "Observed"))) %>%
  select(Status, Observed = x2)

# imp methods

## mean
imps_mean <- mice(
  dat_train,
  m = 1,
  formulas = formulas,
  method = "mean",
  printFlag = FALSE,
  seed = 1,
  maxit = 1
)

# regression
imps_reg <- mice(
  dat_train,
  m = 1,
  formulas = formulas,
  method = "norm.predict",
  printFlag = FALSE,
  seed = 1,
  maxit = 1
)

# stochastic regression
imps_streg <- mice(
  dat_train,
  m = 1,
  formulas = formulas,
  method = "norm.nob",
  printFlag = FALSE,
  seed = 1
)

# boot
imps_boot <- mice(
  dat_train,
  m = 1,
  formulas = formulas,
  method = "norm.boot",
  printFlag = FALSE,
  seed = 1
)

# pmm1-t0
imps_pmm1_t0 <- mice(
  dat_train,
  m = 1,
  formulas = formulas,
  method = "pmm",
  printFlag = FALSE,
  donors = 1,
  seed = 1,
  maxit = 1,
  matchtype = 0
)

# pmm3-t0
imps_pmm3_t0 <- mice(
  dat_train,
  m = 1,
  formulas = formulas,
  method = "pmm",
  printFlag = FALSE,
  donors = 3,
  seed = 1,
  matchtype = 0
)

# pmm1
imps_pmm1 <- mice(
  dat_train,
  m = 1,
  formulas = formulas,
  method = "pmm",
  printFlag = FALSE,
  donors = 1,
  seed = 1
)

# pmm3
imps_pmm3 <- mice(
  dat_train,
  m = 1,
  formulas = formulas,
  method = "pmm",
  printFlag = FALSE,
  donors = 3,
  seed = 1
)

dat_plot <- dat_plot %>%
  mutate(
    x2_mean = complete(imps_mean, 1)$x2,
    x2_reg = complete(imps_reg, 1)$x2,
    x2_streg = complete(imps_streg, 1)$x2,
    x2_boot = complete(imps_boot, 1)$x2,
    x2_pmm1_t0 = complete(imps_pmm1_t0, 1)$x2,
    x2_pmm3_t0 = complete(imps_pmm3_t0, 1)$x2,
    x2_pmm1 = complete(imps_pmm1, 1)$x2,
    x2_pmm3 = complete(imps_pmm3, 1)$x2,
  ) %>%
  mutate(across(c(x2_mean, x2_reg, x2_streg, x2_boot, x2_pmm1_t0, x2_pmm3_t0, x2_pmm1, x2_pmm3), \(x) if_else(Status == "Observed", NA, x))) %>%
  mutate(x2_bayes = NA)

# bayesian
bform <- bf(y | mi() ~ x1 + mi(x2)) +
  bf(x2 | mi() ~ 1) + set_rescor(FALSE)
full_bayes <- brm(bform, data = dat_train, cores = 1)
x2_bayes_miss <- as_draws_df(full_bayes) %>%
  select(contains("Ymi_x2")) %>%
  slice(3909) %>%
  unlist() %>%
  unname()

dat_plot$x2_bayes[dat_plot$Status == "Imputed"] <- x2_bayes_miss

dat_plot_new <- dat_plot %>%
  select(
    Status,
    Observed,
    Mean = x2_mean,
    Reg = x2_reg,
    StReg = x2_streg,
    Boot = x2_boot,
    `PMM-1D-T0` = x2_pmm1_t0,
    `PMM-3D-T0` = x2_pmm3_t0,
    `PMM-1D-T1` = x2_pmm1,
    `PMM-3D-T1` = x2_pmm3,
    FBDA = x2_bayes
  ) %>%
  pivot_longer(cols = c(-Status), names_to = "Type", values_to = "x2", values_drop_na = TRUE) %>%
  mutate(
    Type = factor(Type, levels = rev(c("Observed", "Mean", "Reg", "StReg", "Boot", "PMM-1D-T0", "PMM-3D-T0", "PMM-1D-T1", "PMM-3D-T1", "FBDA"))),
    Status = factor(Status, levels = c("Observed", "Imputed"))
    )


imp_plot <- ggplot(dat_plot_new, aes(x = Type, y = x2, color = Status)) +
  geom_point(alpha = 0.4) +
  scale_color_viridis_d(option = "B", begin = 0.8, end = 0.2) +
  lims(y = c(-3, 3)) +
  labs(y = "x2 Explanatory Variable", x = "Missing Data Approach") +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

imp_plot


# default 6.03 x 4.97 in
ggsave(
  filename = here("inst", "figures", "imputations-0.jpeg"),
  plot = imp_plot,
  dpi = 300,
  units = "in",
  width = 6.03,
  height = 4.97
)
