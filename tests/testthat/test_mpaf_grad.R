library(testthat)

context("Mortality PAF gradients")
library(rpaf)

test_that("hazard and survival gradients", {
  z <- matrix(scan(text ="
                   1 0 0
                   0 1 0
                   1 0 0
                   0 1 0
                   1 0 1
                   0 1 1
                   1 0 1
                   0 1 1", quiet = TRUE), ncol = 3, byrow = TRUE)
  ID <- gl(4, 2)
  PERIOD <- gl(2, 1, 8)
  breaks = 0:2
  lambda <- c(2,3,2,3,6,7,6,7)

  S <- exp(-lambda)

  # grad <- dpaf_grad(z, hz_d, hz_m, sv_d, sv_m, breaks, ID, PERIOD)

  # expected values
  e_glambda <- matrix(c(2,0,0, 0,3,0, 2,0,0, 0,3,0, 6,0,6, 0,7,7, 6,0,6, 0,7,7),
                   ncol = 3, byrow = TRUE)
  e_glambda_sum <- matrix(c(2,0,0, 2,3,0, 2,0,0, 2,3,0,
                        6,0,6, 6,7,13, 6,0,6, 6,7,13),
                      ncol = 3, byrow = TRUE)
  e_gS <- -exp(-lambda) * e_glambda_sum

  glambda <- rpaf:::mpaf_grad_lambda(z, lambda)
  expect_equal(glambda, e_glambda)

  gS <- rpaf:::mpaf_grad_S(glambda, S, ID, PERIOD, diff(breaks))
  expect_equal(gS, e_gS)
})

test_that("I gradients, single-person", {
  z <- diag(3) # 3x3 identity matrix
  id <- rep(1, 3)
  period <- gl(3, 1)
  dt <- c(1,1,1)

  S <- c(0.99, 0.95, 0.90)
  lambda <- log(c(1/0.99, 0.99/0.95, 0.95/0.90))

  stopifnot(all.equal(S, rpaf:::mpaf_S(lambda, id, period, dt)))

  glambda <- diag(lambda)

  stopifnot(all.equal(glambda, rpaf:::mpaf_grad_lambda(z, lambda)))

  gS <- -S %o% lambda
  gS[upper.tri(gS)] <- 0

  expect_equal(gS, rpaf:::mpaf_grad_S(glambda, S, id, period, dt))

  gdS <- apply(gS, 2, function(col) -diff(c(0, col)))

  expect_equal(gdS, rpaf:::mpaf_grad_I(gdS, id, period),
               check.attributes = FALSE)
})

test_that("I gradients, single-period", {
  z <- diag(3) # 3x3 identity matrix
  id <- rep(3, 1)
  period <- gl(1, 3)
  dt <- c(1)

  S <- c(0.99, 0.95, 0.90)
  lambda <- log(c(1/0.99, 1/0.95, 1/0.90))

  stopifnot(all.equal(S, rpaf:::mpaf_S(lambda, id, period, dt)))

  glambda <- diag(lambda)

  stopifnot(all.equal(glambda, rpaf:::mpaf_grad_lambda(z, lambda)))

  gS <- -S %o% lambda
  gS[upper.tri(gS)] <- 0

  expect_equal(gS, rpaf:::mpaf_grad_S(glambda, S, id, period, dt))

  gdS <- apply(gS, 2, function(col) -diff(c(0, col)))

  expect_equal(gdS, rpaf:::mpaf_grad_I(gdS, id, period),
               check.attributes = FALSE)
})

