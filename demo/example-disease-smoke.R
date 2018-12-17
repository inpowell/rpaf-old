# Example user script for calculating PAF for disease
#
# This script only refers to functions in R/disease_user.R -- see that source
# file, and the documentation, for more information.
#
# In terms of SAS examples, this aims to replicate the `paf_example_d_smoke`
# output.

library(rpaf)

# prepare the data, similar in function to the SAS gen_data macro, but
# specialised for disease PAFs
smoke_data <- dpaf_data(
  rawdata = minifhs,
  variables = c("B_COHORT", "SEX", "SMOKE"),
  time_breaks = c(0,5,10,15,17), # could use, e.g. `c(seq(0,17,5), 17)`
  id_var = "ID",
  death_time = "DEATH_FT",
  death_ind = "DEATH",
  disease_time = "DIAB_FT",
  disease_ind = "DIAB",
  period_factor = "F_PERIOD",
  ft_end = "F_TIME"
)

# Fit the survival regression models, and return the information obatined. This
# is similar in function to the first part of the est_matrix macro, but doesn't
# calculate relative risks. Instead, it only fits the survival regression models
# for both disease and mortality, and returns the model frame (which can be used
# to prepare the design matrix).
smoke_fit <- dpaf(
  disease_resp = survival::Surv(F_TIME, DIAB) ~ .,
  death_resp = survival::Surv(F_TIME, DEATH) ~ .,
  predictors = ~ 0 + B_COHORT*F_PERIOD + SEX + SMOKE,
  dpaf_data = smoke_data
)

# internally, this calls the summary.dpaf function in R/disease_user.R, but is
# wrapped around a well-known R function.
smoke_summ <- summary(smoke_fit,
                      modlist = list(SMOKE = c("Never", "<30/day", ">=30/day")),
                      comprehensive = TRUE)

smoke_summ
