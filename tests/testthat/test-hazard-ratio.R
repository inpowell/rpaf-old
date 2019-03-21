context("test-hazard-ratio")

library(testthat)

df <- data.frame(
  tmt = gl(2, 1, labels = c("N", "Y")),
  died = c(TRUE, TRUE),
  time = c(1.0, 2.0),
  id = 1:2
)

test_that("hazard ratio reported correctly", {
  #  expect_true(FALSE && "TODO")
  sr <- survival::survreg(survival::Surv(time, died) ~ tmt + 1, data = df, dist = "exp")
  hr <- rpaf::hazard_ratios(sr, "tmt")
  expect_equal(
    hr$tmt[[1]],
    rbind(c(1,1,1), exp(log(0.5) + sqrt(2) * qnorm(c(0.5, 0.025, 0.975)))),
    check.attributes = FALSE
  )
  expect_equal(names(hr), c("tmt"))
  expect_equal(dimnames(hr$tmt), list(ALL = ""))
})
