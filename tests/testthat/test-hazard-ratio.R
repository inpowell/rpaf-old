context("test-hazard-ratio")

library(testthat)

easy_dat <- data.frame(
  tmt = gl(2, 1, labels = c("N", "Y")),
  died = c(TRUE, TRUE),
  time = c(1.0, 2.0),
  id = 1:2
)

# interaction examples
int_dat <- data.frame(
  tmt = gl(2, 1, 8, labels = c("N", "Y")),
  blk = gl(2, 2, 8, labels = c("A", "B")),
  died = rep(TRUE, 4),
  time = c(1.0, 2.0, 3.0, 4.0, 2.0, 3.0, 4.0, 5.0),
  id = 1:8
)
int_sr <- survival::survreg(survival::Surv(time, died) ~ tmt * blk + 1,
                            data = int_dat, dist = "exp")
int_hr <- rpaf::hazard_ratios(int_sr, "tmt")

test_that("hazard ratio reported correctly", {
  #  expect_true(FALSE && "TODO")
  sr <- survival::survreg(survival::Surv(time, died) ~ tmt + 1,
                          data = easy_dat, dist = "exp")
  hr <- rpaf::hazard_ratios(sr, "tmt")
  expect_equal(
    hr$tmt[[1]],
    rbind(c(1,1,1), exp(log(0.5) + sqrt(2) * qnorm(c(0.5, 0.025, 0.975)))),
    check.attributes = FALSE
  )
  expect_equal(names(hr), c("tmt"))
  expect_equal(dimnames(hr$tmt), list(ALL = ""))
})

test_that("interactions give correct row names of hazard ratios", {
  # expect_equal(names(hr), c("tmt"))
  expect_equal(dimnames(int_hr$tmt), list(blk = c("A", "B")))
  expect_equal(rownames(int_hr$tmt[[1]]), c("N", "Y"))
})

test_that("first row of hazard ratio tables is always 1", {
  expect_equivalent(int_hr$tmt[[1]][1,], c(1.0, 1.0, 1.0))
  expect_equivalent(int_hr$tmt[[2]][1,], c(1.0, 1.0, 1.0))
})
