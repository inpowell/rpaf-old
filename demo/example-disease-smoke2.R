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
smoke_data <- gen_data(
  indata = minifhs,
  id_var = "ID",
  ft_breaks = c(0,5,10,15,17),
  death_ind = "DEATH",
  death_time = "DEATH_FT",
  disease_ind = "DIAB",
  disease_time = "DIAB_FT",
  variables = c("B_COHORT", "SEX", "SMOKE"),
  period_factor = "F_PERIOD",
  time_var = "F_TIME"
)

smoke_fitd <- est_matrix(survival::Surv(F_TIME, DIAB) ~
                           B_COHORT * F_PERIOD + SEX + SMOKE,
                         smoke_data,
                         list(SMOKE = c("Never", "<30/day", ">=30/day")),
                         covar_model = c("SEX", "SMOKE"))
smoke_fitm <- est_matrix(survival::Surv(F_TIME, DEATH) ~
                           B_COHORT * F_PERIOD + SEX + SMOKE,
                         smoke_data,
                         list(SMOKE = c("Never", "<30/day", ">=30/day")),
                         covar_model = c("SEX", "SMOKE"))

smoke_paf <- dpaf_est_paf(smoke_fitd, smoke_fitm, smoke_data)

smoke_summ <- dpaf_summary(
  disease_resp = survival::Surv(F_TIME, DIAB) ~ .,
  death_resp = survival::Surv(F_TIME, DEATH) ~ .,
  predictors = ~ 0 + B_COHORT * F_PERIOD + SEX + SMOKE,
  dpaf_data = smoke_data,
  modifications = list(SMOKE = c("Never", "<30/day", ">=30/day")),
  covar_model = c("SEX", "SMOKE")
)
