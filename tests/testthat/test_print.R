library(testthat)

context("Printing functions")

library(rpaf)

df <- minifhs[c(1:20, 1096),]
smoke_data <- gen_data(df, id_var = "ID", ft_breaks = c(0,17),
                       primary_ind = "DEATH", primary_time = "DEATH_FT",
                       variables = c("B_COHORT", "SEX", "SMOKE"))

smoke_mpaf <- mpaf_summary(
  survival::Surv(f_end, DEATH) ~ B_COHORT + SEX + SMOKE,
  smoke_data,
  modifications = list(SMOKE = c("Never", "<30/day", ">=30/day")),
  covar_model = c("SEX", "SMOKE")
)

bmi_data <- gen_data(df, id_var = "ID", ft_breaks = c(0,17),
                     primary_ind = "DIAB", primary_time = "DIAB_FT",
                     secondary_ind = "DEATH", secondary_time = "DEATH_FT",
                     variables = c("B_COHORT", "SEX", "BMI_2")
                    )
bmi_dpaf <- dpaf_summary(
  survival::Surv(f_end, DIAB) ~ .,
  survival::Surv(f_end, DEATH) ~ .,
  ~ B_COHORT + SEX + BMI_2,
  bmi_data,
  modifications = list(BMI_2 = "<25.0"),
  covar_model = c("SEX", "BMI_2")
)

test_that("disease PAF summary prints output headings", {
  expect_output(print(bmi_dpaf), "Call:\ndpaf_summary", fixed = TRUE)
  # expect_output(print(bmi_dpaf), "Data preparation call:\ngen_data",
  #               fixed = TRUE)
  expect_output(print(bmi_dpaf), "Survival regression summary:", fixed = TRUE)
  expect_output(print(bmi_dpaf),
                "Cohort size: 20 (1 observation deleted due to missingness)",
                fixed = TRUE)
  expect_output(print(bmi_dpaf), "Hazard ratios:", fixed = TRUE)
  expect_output(print(bmi_dpaf), "Modifications:\nBMI_2: . -> <25.0")
  expect_output(print(bmi_dpaf), "PAFs for disease:")
})

test_that("mortality PAF summary prints output headings", {
  expect_output(print(smoke_mpaf), "Call:\nmpaf_summary", fixed = TRUE)
  # expect_output(print(smoke_mpaf), "Data preparation call:\ngen_data",
  #               fixed = TRUE)
  expect_output(print(smoke_mpaf), "Survival regression summary:", fixed = TRUE)
  expect_output(print(smoke_mpaf), "Cohort size: 21", fixed = TRUE)
  expect_output(print(smoke_mpaf), "Hazard ratios:", fixed = TRUE)
  expect_output(print(smoke_mpaf), "Modifications:\nSMOKE: <30/day, >=30/day -> Never")
  expect_output(print(smoke_mpaf), "PAFs for death:")
})
