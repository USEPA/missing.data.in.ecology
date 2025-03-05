library(tidyverse)
library(here)
library(patchwork)
set.seed(1)
n <- 20000

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


######################## MCAR
p <- 0.5
VMMI <- rnorm(2 * n, mean = 38, sd = 14)
miss <- sample(c("Missing", "Observed"), size = 2 * n, prob = c(p, 1 - p), replace = TRUE)
miss <- factor(miss, levels = c("Observed", "Missing"))
mcar_dat <- tibble(VMMI = VMMI, miss = miss, type = "All")

p_mcar <- ggplot(mcar_dat, aes(x = VMMI)) +
    geom_density(color = okabe[1], bw = 4, linewidth = 0.8) +
    facet_grid(~ miss) +
    labs(x = "VMMI", y = "Density", title = "a)") +
    lims(x = c(0, 100)) +
    theme_bw() +
    theme(
      axis.text.y = element_blank()
    )


p_mcar

######################## MAR
p <- 0.7
wade_VMMI <- rnorm(n, mean = 50, sd = 12)
wade_miss <- sample(c("Missing", "Observed"), size = n, prob = c(p, 1 - p), replace = TRUE)
boat_VMMI <- rnorm(n, mean = 30, sd = 8)
boat_miss <- sample(c("Missing", "Observed"), size = n, prob = c(1 - p, p), replace = TRUE)
VMMI <- c(wade_VMMI, boat_VMMI)
miss <- c(wade_miss, boat_miss)
miss <- factor(miss, levels = c("Observed", "Missing"))
type <- rep(c("Herbaceous", "Woody"), each = n)
mar_dat <- tibble(type = type, VMMI = VMMI, miss = miss)
mar_dat <- bind_rows(mar_dat, mar_dat %>% mutate(type = "All"))

p_mar <- ggplot(mar_dat, aes(x = VMMI, color = type, linetype = type)) +
  geom_density(bw = 3, linewidth = 0.8) +
  facet_grid(~ miss) +
  scale_color_manual(values = okabe[1:3], name = "Wetland Type") +
  scale_linetype_manual(values = c(1, 2, 5), name = "Wetland Type") +
  labs(x = "VMMI", y = "Density", title = "b)") +
  lims(x = c(0, 100)) +
  theme_bw() +
  theme(
    axis.text.y = element_blank()
  )

######################## MNAR
soilhard1 <- c(runif(n/2, 0, 0), runif(n/2, 1, 1))
wade_VMMI <- c(rnorm(n/2, mean = 34, sd = 13), rnorm(n/2, mean = 68, sd = 7))
wade_miss <- rbinom(n, size = 1, prob = soilhard1)
wade_miss <- ifelse(wade_miss == 0, "Missing", "Observed")
soilhard2 <- c(runif(n/2, 0, 0), runif(n/2, 1, 1))
boat_VMMI <- c(rnorm(n/2, mean = 22, sd = 6), rnorm(n/2, mean = 68, sd = 14))
boat_miss <- rbinom(n, size = 1, prob = soilhard2)
boat_miss <- ifelse(boat_miss == 0, "Missing", "Observed")
VMMI <- c(wade_VMMI, boat_VMMI)
miss <- c(wade_miss, boat_miss)
miss <- factor(miss, levels = c("Observed", "Missing"))
type <- rep(c("Herbaceous", "Woody"), each = n)
nmar_dat <- tibble(type = type, VMMI = VMMI, miss = miss)
nmar_dat <- bind_rows(nmar_dat, nmar_dat %>% mutate(type = "All"))

p_nmar <- ggplot(nmar_dat, aes(x = VMMI, color = type, linetype = type)) +
  geom_density(bw = 3, linewidth = 0.8) +
  facet_grid(~ miss) +
  scale_color_manual(values = okabe[1:3], name = "Wetland Type") +
  scale_linetype_manual(values = c(1, 2, 5), name = "Wetland Type") +
  labs(x = "VMMI", y = "Density", title = "c)") +
  lims(x = c(0, 100)) +
  theme_bw() +
  theme(
    axis.text.y = element_blank()
  )

p1 <- p_mcar + p_mar + p_nmar + plot_layout(nrow = 3, guides = "collect") # & theme(legend.position = "top")
p1

ggsave(
  filename = here("inst", "figures", "patterns.jpeg"),
  plot = p1,
  dpi = 300,
  units = "in"
)

