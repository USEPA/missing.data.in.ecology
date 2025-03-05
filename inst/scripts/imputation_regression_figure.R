library(tidyverse)
library(mice)
library(patchwork)
library(here)


# use patchwork
set.seed(57) #25, 57
x1 <- seq(0, 1, length.out = 41)
x2 <- 0 + 3 * x1 + rnorm(length(x1), sd = 1)
imp <- ifelse(as.character(x1) %in% as.character(c(seq(0.05, 0.95, by = 0.1))), "yes", "no")
x2_imp <- ifelse(imp == "yes", NA, x2)
size <- ifelse(imp == "yes", 3, 1.5)
dat <- data.frame(x1, x2, imp, x2_imp)
formulas <- list(x2_imp ~ x1)
formulas <- name.formulas(formulas)

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# panel A
p1 <- ggplot(dat, aes(x = x1, y = x2_imp)) +
  geom_smooth(method = "lm", se = FALSE, color = okabe[3]) +
  geom_point(color = okabe[2], size = 1.5) +
  labs(x = expression(z[1]), y = expression(z[2]), title = "Observed Data") +
  lims(y = c(-2, 4.2)) +
  scale_y_continuous(limits = c(-2, NA), breaks = seq(-2, 4, length.out = 4)) +
  theme_bw(base_size = 10) +
  theme(
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text=element_text(size=14)
  )




# panel B
imps <- mice(
  dat,
  m = 1,
  formulas = formulas,
  method = "mean",
  printFlag = FALSE
)
dat_imp_sub <- complete(imps, 1)
p2 <- ggplot(dat_imp_sub, aes(x = x1, y = x2_imp, color = imp, shape = imp)) +
  geom_smooth(
    data = dat_imp_sub %>% filter(imp == "no"),
    method = "lm",
    se = FALSE,
    color = okabe[3]
  ) +
  geom_point(size = size) +
  scale_color_manual(name = "Imputed", values = okabe[1:2], breaks = c("yes", "no"), labels = c("yes", "no")) +
  scale_shape_manual(name = "Imputed", values = c(17, 19), breaks = c("yes", "no"), labels = c("yes", "no")) +
  scale_y_continuous(limits = c(-2, NA), breaks = seq(-2, 4, length.out = 4)) +
  labs(x = expression(z[1]), y = expression(z[2]), title = "Mean Imputation") +
  theme_bw(base_size = 10) +
  theme(
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text=element_text(size=12)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2)))


# panel C
imps <- mice(
  dat,
  m = 1,
  formulas = formulas,
  method = "norm.predict",
  printFlag = FALSE,
  maxit = 1
)
dat_imp_sub <- complete(imps, 1)
p3 <- ggplot(dat_imp_sub, aes(x = x1, y = x2_imp, color = imp, shape = imp)) +
  geom_smooth(
    data = dat_imp_sub %>% filter(imp == "no"),
    method = "lm",
    se = FALSE,
    color = okabe[3]
  ) +
  geom_point(size = size) +
  scale_color_manual(name = "Imputed", values = okabe[1:2], breaks = c("yes", "no"), labels = c("yes", "no")) +
  scale_shape_manual(name = "Imputed", values = c(17, 19), breaks = c("yes", "no"), labels = c("yes", "no")) +
  scale_y_continuous(limits = c(-2, NA), breaks = seq(-2, 4, length.out = 4)) +
  labs(x = expression(z[1]), y = expression(z[2]), title = "Reg Imputation") +
  theme_bw(base_size = 10) +
  theme(
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text=element_text(size = 12)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

# panel D
p4 <- imps <- mice(
  dat,
  m = 1,
  formulas = formulas,
  method = "norm.boot",
  printFlag = FALSE,
  maxit = 5
)
dat_imp_sub <- complete(imps, 1)
p4 <- ggplot(dat_imp_sub, aes(x = x1, y = x2_imp, color = imp, shape = imp)) +
  geom_smooth(
    data = dat_imp_sub %>% filter(imp == "no"),
    method = "lm",
    se = FALSE,
    color = okabe[3]
  ) +
  geom_point(size = size) +
  scale_color_manual(name = "Imputed", values = okabe[1:2], breaks = c("yes", "no"), labels = c("yes", "no")) +
  scale_shape_manual(name = "Imputed", values = c(17, 19), breaks = c("yes", "no"), labels = c("yes", "no")) +
  scale_y_continuous(limits = c(-2, NA), breaks = seq(-2, 4, length.out = 4)) +
  labs(x = expression(z[1]), y = expression(z[2]), title = "Boot Imputation") +
  theme_bw(base_size = 10) +
  theme(
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text=element_text(size=12)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

p <- p1 + p2 + p3 + p4 + plot_layout(guides = "collect")
p
ggsave(
  filename = here("inst", "figures", "imputations-comp-0.jpeg"),
  plot = p,
  dpi = 300,
  height = 4.98,
  width = 6.39
)
