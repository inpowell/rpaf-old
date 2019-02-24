library(testthat)

context("Model fitting and hazard ratios")

library(rpaf)

df <- data.frame(
  tmt = gl(2, 1, labels = c("N", "Y")),
  died = c(TRUE, TRUE),
  time = c(1.0, 2.0),
  id = 1:2
)

paf_data <- gen_data(df, id_var = "id", ft_breaks = c(0,5),
                     primary_ind = "died", primary_time = "time", variables = "tmt",
                     time_var = "f_time"
                     )
test_that("reference levels", {
  paf_fit <- est_matrix(
    survival::Surv(f_time, died) ~ 1 + tmt, paf_data,
    modifications = list(tmt = "Y"),
    covar_model = c("tmt")
  )
  expect_equal(paf_fit$HR[,1], c("tmtN" = 1, "tmtY" = 0.5))
  expect_equal(is.na(paf_fit$HR[1,]), c(FALSE, TRUE, TRUE),
               check.attributes = FALSE)
})
