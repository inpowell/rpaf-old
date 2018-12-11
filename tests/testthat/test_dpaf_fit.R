library(testthat)
context("Morbidity")

library(rpaf)

test_that("Single-period", {
  x <- matrix(c(1,0,0, 0,1,1), ncol = 2)
  cf_d <- c(g1 = -log(0.5), g2 = -log(3))
  cf_m <- c(g1 = -log(2), g2 = -1)

  id <- gl(3,1)
  period <- gl(1,3, ordered = TRUE)

  ss <- dpaf_fit(x, cf_d, cf_m, id, period, 0:1)

  expect_equal(ss$hazard.disease, c(0.5, 3, 3))
  expect_equal(ss$hazard.mortality, c(2, exp(1), exp(1)))
  expect_equal(ss$survival.disease, c(exp(-0.5), exp(-3), exp(-3)))
  expect_equal(ss$survival.mortality, c(exp(-2), exp(-exp(1)), exp(-exp(1))))
  expect_equal(
    ss$morbidity,
    ( (1 - exp(-2.5)) / 5 +
        2 * (3 / (3 + exp(1))) * (1 - exp(-3-exp(1))) ) / 3
  )
})

test_that("Single-person", {
  x <- matrix(c(1,0,0, 0,1,1), ncol = 2)
  cf_d <- c(g1 = -log(0.5), g2 = -log(3))
  cf_m <- c(g1 = -log(2), g2 = -log(5))

  id <- gl(1,3)
  period <- gl(3,1, ordered = TRUE)

  ss <- dpaf_fit(x, cf_d, cf_m, id, period, 0:3)

  expect_equal(ss$hazard.disease, c(0.5, 3, 3))
  expect_equal(ss$hazard.mortality, c(2, 5, 5))
  expect_equal(ss$survival.disease, c(exp(-0.5), exp(-3.5), exp(-6.5)))
  expect_equal(ss$survival.mortality, c(exp(-2), exp(-7), exp(-12)))
  expect_equal(
    ss$morbidity,
    (1 - exp(-2.5))/5 +
      (exp(-2.5) - exp(-10.5)) * 3 / 8 +
      (exp(-10.5) - exp(-18.5)) * 3 / 8
  )
})

test_that("Uneven period", {
  x <- matrix(c(1,0,0, 0,1,1), ncol = 2)
  cf_d <- c(g1 = -log(0.5), g2 = -log(3))
  cf_m <- c(g1 = -log(2), g2 = -log(3))

  id <- gl(1, 3)
  period <- gl(3, 1, ordered = TRUE)

  ss <- dpaf_fit(x, cf_d, cf_m, id, period, c(0, 1, 2, 7/3))

  expect_equal(ss$hazard.disease, c(0.5, 3, 3))
  expect_equal(ss$hazard.mortality, c(2, 3, 3))
})

test_that("Unsorted periods", {
  x <- matrix(c(1,0,0, 0,1,1), ncol = 2)
  cf_d <- c(g1 = -log(0.5), g2 = -log(3))
  cf_m <- c(g1 = -log(2), g2 = -log(3))

  id <- gl(1, 3)
  period <- factor(c(1,3,2), levels = 1:3, ordered = TRUE)
  period2 <- sort(period)
  x2 <- x[order(period),]
  id2 <- id[order(period)]

  expect_warning(f <- dpaf_fit(x, cf_d, cf_m, id, period, 0:3),
                 "not correctly sorted")
  expect_equal(f, dpaf_fit(x2, cf_d, cf_m, id2, period2, 0:3))
})
