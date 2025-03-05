bold_y <- function() {
  y_1 <- rnorm(n = 6, mean = 0, sd = 1)
  y_2 <- rpois(n = 6, lambda = 2)
  y <- matrix(c(y_1, y_2), nrow = 6, ncol = 2)
  return(y)
}
approx_bold_Y <- function(s) {
  Y <- replicate(s, bold_y())
  return(Y)
}

set.seed(1)
y_draw_one <- bold_y()
print(y_draw_one)
y_draw_two <- bold_y()
print(y_draw_two)

set.seed(1)
approx_bold_Y_draw <- approx_bold_Y(s = 10000)
print(approx_bold_Y_draw[, , 1])
print(approx_bold_Y_draw[, , 10000])
head(approx_bold_Y_draw[1, 1, ], n = 40)


bold_r <- function() {
  r_1 <- rbinom(n = 6, size = 1, prob = 0.9)
  r_2 <- rbinom(n = 6, size = 1, prob = 0.7)
  r <- matrix(c(r_1, r_2), nrow = 6, ncol = 2)
  return(r)
}
approx_bold_R <- function(s) {
  R <- replicate(s, bold_r())
  return(R)
}
set.seed(1)
r_draw_one <- bold_r()
print(r_draw_one)
r_draw_two <- bold_r()
print(r_draw_two)

set.seed(1)
approx_bold_R_draw <- approx_bold_R(s = 10000)
print(approx_bold_R_draw[, , 1])
print(approx_bold_R_draw[, , 10000])
head(approx_bold_R_draw[1, 1, ], n = 40)
