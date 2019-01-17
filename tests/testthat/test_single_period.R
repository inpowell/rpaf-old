library(testthat)

context("Single-period PAF studies")

library(rpaf)

gdata <- gen_data(
  indata = head(minifhs, 20), # contains deaths and one diabetes case
  id_var = "ID", ft_breaks = c(0,17),
  death_ind = "DEATH", death_time = "DEATH_FT",
  disease_ind = "DIAB", disease_time = "DIAB_FT",
  variables = c("B_COHORT", "SEX", "BMI_2")
) # strictly it's not okay to use this for both disease and death

test_that("single-period mortality", {
  summary_m <- mpaf_summary(
    survival::Surv(f_end, DEATH) ~ B_COHORT + SEX + BMI_2,
    mpaf_data = gdata, modifications = list(BMI_2 = "<25.0"),
    covar_model = c("SEX", "BMI_2")
  )

  expect_equal(dimnames(summary_m$paf0),
               list("(0,17]", c("PAF", "2.5 %", "97.5 %")))
})

test_that("single-period disease", {
  summary_d <- dpaf_summary(
    survival::Surv(f_end, DIAB) ~ .,
    survival::Surv(f_end, DEATH) ~ .,
    ~ B_COHORT + SEX + BMI_2,
    dpaf_data = gdata, modifications = list(BMI_2 = "<25.0"),
    covar_model = c("SEX", "BMI_2")
  )

  expect_equal(dimnames(summary_d$paf0),
               list("(0,17]", c("PAF", "2.5 %", "97.5 %")))
})
