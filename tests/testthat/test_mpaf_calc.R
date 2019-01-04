library(testthat)

context("Mortality PAF calculations")

library(rpaf)

test_that("Single-period", {
  x <- matrix(c(1,0,0, 0,1,1), ncol = 2)
  cf <- c(g1 = -log(2), g2 = -1)

  id <- gl(3,1)
  period <- gl(1,3, ordered = TRUE)

  lda <- rpaf:::mpaf_lambda(x, cf)
  S <- rpaf:::mpaf_S(lda, id, period, dt = c(1))
  dS <- rpaf:::mpaf_delta_S(S, id, period)
  I <- rpaf:::mpaf_I(dS, period)

  expect_equal(lda, c(2, exp(1), exp(1)))
  expect_equal(S, c(exp(-2), exp(-exp(1)), exp(-exp(1))))
  expect_equal(dS, 1 - S)
  expect_equal(I, 1 - (exp(-2) + 2*exp(-exp(1))) / 3)
})

test_that("Single-person", {
  x <- matrix(c(1,0,0, 0,1,1), ncol = 2)
  cf <- c(g1 = -log(2), g2 = -log(5))

  id <- gl(1,3)
  period <- gl(3,1)

  lambda <- rpaf:::mpaf_lambda(x, cf)
  S <- rpaf:::mpaf_S(lambda, id, period, c(1,1,1))
  dS <- rpaf:::mpaf_delta_S(S, id, period)
  I <- rpaf:::mpaf_I(dS, period)

  expect_equal(lambda, c(2, 5, 5))
  expect_equal(S, c(exp(-2), exp(-7), exp(-12)))
  expect_equal(dS, c(1 - exp(-2), exp(-2) - exp(-7), exp(-7) - exp(-12)))
  expect_equal(I, c(1 - exp(-2), exp(-2) - exp(-7), exp(-7) - exp(-12)))
})

