library(mice)
set.seed(1)
# regression imputation

vmmi <- c(44, 31.2, 39.4, 25.0, 41.2, 38.7)
vmmi_m <- vmmi
miss_index <- c(2, 6)
vmmi_m[miss_index] <- NA
mean(vmmi_m, na.rm = TRUE)
type <- c("W", "H", "H", "W", "H", "W")
hard <- c(97.5, 29.5, 90.4, 79.6, 33.6, 65.9)

dat <- data.frame(vmmi, vmmi_m, type, hard)

# linear model
lmod <- lm(vmmi_m ~ type + hard, data = dat)
summary(lmod)
coef(lmod)
preds <- predict(lmod, newdata = dat[miss_index, ])
preds

formulas <- list(
  vmmi_m ~ type + hard
)
formulas <- name.formulas(formulas)
# matches above
imps <- mice(dat[, -1], formulas = formulas, method = "mean", m = 1, maxit = 1)
complete(imps, 1)

imps <- mice(dat[, -1], formulas = formulas, method = "norm.predict", m = 1, maxit = 1)
complete(imps, 1)

errs <- rnorm(n = 2, sd = summary(lmod)$sigma)
errs
preds + errs
summary(lmod)$sigma

# bootstrap
lmod <- lm(vmmi_m ~ type + hard, data = dat[c(3, 3, 4, 5), ])
coef(lmod)
preds <- predict(lmod, dat[miss_index, ])
errs <- c(5.14, -2.76)
preds + errs
