library(testthat)
context("Morbidity")

library(rpaf)

test_that("Single-period", {
  x <- matrix(c(1,0,0, 0,1,1), ncol = 2)
  cf_d <- c(g1 = -log(0.5), g2 = -log(3))
  cf_m <- c(g1 = -log(2), g2 = -1)

  id <- gl(3,1)
  period <- gl(1,3, ordered = TRUE)

  hz <- rpaf:::dpaf_hz(x, cbind("disease" = cf_d, "mortality" = cf_m))
  sv <- rpaf:::dpaf_sv(hz, id, period, dt = c(1))

  expect_equal(hz, cbind(disease = c(0.5, 3, 3),
                         mortality = c(2, exp(1), exp(1))))
  expect_equal(sv, cbind(disease = c(exp(-0.5), exp(-3), exp(-3)),
                         mortality = c(exp(-2), exp(-exp(1)), exp(-exp(1)))))
  expect_equal(
    rpaf:::dpaf_i(hz, rpaf:::dpaf_svp(sv), id, period),
    0.5 / 2.5 * (1 - exp(-2.5)) + 3 / (3 + exp(1)) * (1 - exp(-3-exp(1))) * 2
  )
})

test_that("Single-person", {
  x <- matrix(c(1,0,0, 0,1,1), ncol = 2)
  cf_d <- c(g1 = -log(0.5), g2 = -log(3))
  cf_m <- c(g1 = -log(2), g2 = -log(5))
  cf <- cbind(disease = cf_d, mortality = cf_m)

  id <- gl(1,3)
  period <- gl(3,1)

  # ss <- dpaf_fit(x, cf_d, cf_m, id, period, 0:3)

  hz <- rpaf:::dpaf_hz(x, cf)
  sv <- rpaf:::dpaf_sv(hz, id, period, c(1,1,1))

  expect_equal(hz, cbind(disease = c(0.5, 3, 3),
                         mortality = c(2, 5, 5)))
  expect_equal(sv, cbind(disease = c(exp(-0.5), exp(-3.5), exp(-6.5)),
                         mortality = c(exp(-2), exp(-7), exp(-12))))

  expect_equal(
    rpaf:::dpaf_i(hz, rpaf:::dpaf_svp(sv), id, period),
    0.5 / 2.5 * (1 - exp(-2.5)) +
      3 / 8 * (exp(-2.5) - exp(-10.5)) +
      3 / 8 * (exp(-10.5) - exp(-18.5))
  )
})

test_that("Uneven period", {
  x <- matrix(c(1,0,0, 0,1,1), ncol = 2)
  cf_d <- c(g1 = -log(0.5), g2 = -log(3))
  cf_m <- c(g1 = -log(2), g2 = -log(5))
  cf <- cbind(disease = cf_d, mortality = cf_m)

  id <- gl(1,3)
  period <- gl(3,1)

  # ss <- dpaf_fit(x, cf_d, cf_m, id, period, 0:3)

  hz <- rpaf:::dpaf_hz(x, cf)
  sv <- rpaf:::dpaf_sv(hz, id, period, c(1,1,1/3))

  expect_equal(hz, cbind(disease = c(0.5, 3, 3),
                         mortality = c(2, 5, 5)))
  expect_equal(sv, cbind(disease = c(exp(-0.5), exp(-3.5), exp(-3.5 - 3/3)),
                         mortality = c(exp(-2), exp(-7), exp(-7 - 5/3))))
})

test_that("Unsorted periods", {
  x <- matrix(c(1,0,0, 0,1,1), ncol = 2)
  cf_d <- c(g1 = -log(0.5), g2 = -log(3))
  cf_m <- c(g1 = -log(2), g2 = -log(3))
  cf <- cbind(disease = cf_d, mortality = cf_m)

  id <- gl(1, 3)
  period <- factor(c(1,3,2), levels = 1:3, ordered = TRUE)

  hz <- rpaf:::dpaf_hz(x, cf)
  expect_error(rpaf:::dpaf_sv(hz, id, period, c(1,1,1/3)))
})
