library(tidyverse)
library(here)
library(spmodel)

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
lmod <- splm(VMMI_2016 ~ PALT_SOHARD + IS_WOODY  + SAMPLEABLE_WATER + SAMPLEABLE_WATER:NTL_RESULT ,
             data = dat, spcov_type = "none")
summary(lmod)
spmod <- splm(VMMI_2016 ~ PALT_SOHARD + IS_WOODY + SAMPLEABLE_WATER + SAMPLEABLE_WATER:NTL_RESULT ,
             data = dat, spcov_type = "exponential", xcoord = "XCOORD", ycoord = "YCOORD")
summary(spmod)
glances(lmod, spmod)
loocv(lmod)
loocv(spmod)
write_csv(tidy(lmod, conf.int = TRUE), str_c(here("inst", "output", "vmmi"), "/fixed_contingency.csv"))
write_csv(tidy(spmod, conf.int = TRUE), str_c(here("inst", "output", "vmmi"), "/fixed_contingency_spmod.csv"))
