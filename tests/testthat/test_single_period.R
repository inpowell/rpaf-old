library(testthat)

context("Single-period PAF studies")

library(rpaf)

gdata <- gen_data(
  indata = head(minifhs, 20), # contains deaths and one diabetes case
  id_var = "ID", ft_breaks = c(0,17),
  primary_ind = "DIAB", primary_time = "DIAB_FT",
  secondary_ind = "DEATH", secondary_time = "DEATH_FT",
  variables = c("B_COHORT", "SEX", "BMI_2")
) # strictly it's not okay to use this for both disease and death

test_that("single-period mortality dimensionality", {
  summary_2 <- mpaf_summary(
    survival::Surv(f_end, DEATH) ~ B_COHORT + SEX + BMI_2,
    mpaf_data = gdata, modifications = list(BMI_2 = "<25.0"),
    covar_model = c("SEX", "BMI_2")
  )

  expect_equal(dimnames(summary_2$paf0),
               list("(0,17]", c("PAF", "2.5 %", "97.5 %")))
})

test_that("single-period disease dimensionality", {
  summary_1 <- dpaf_summary(
    survival::Surv(f_end, DIAB) ~ .,
    survival::Surv(f_end, DEATH) ~ .,
    ~ B_COHORT + SEX + BMI_2,
    dpaf_data = gdata, modifications = list(BMI_2 = "<25.0"),
    covar_model = c("SEX", "BMI_2")
  )

  expect_equal(dimnames(summary_1$paf0),
               list("(0,17]", c("PAF", "2.5 %", "97.5 %")))
})
