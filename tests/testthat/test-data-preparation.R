context("test-data-preparation")

library(testthat)

test_that("collapsing times from multiple events works", {
  testdat <- data.frame(
    ind1 = rep(c(FALSE, TRUE), each = 2, length.out = 8),
    time1 = rep(c(1, 2), each = 4),
    ind2 = rep(c(FALSE, TRUE), each = 1, length.out = 8),
    time2 = rep(c(2, 1), each = 4)
  )

  expdat <- data.frame(
    time = rep(1, 8),
    ind1 = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
    ind2 = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE)
  )

  expect_equal(rpaf::collapse_times(testdat, "ind1", "time1", "ind2", "time2"),
               expdat)
})

test_that("expanding event factor for multiple events works", {
  testdat <- data.frame(
    time = 1:12,
    event = gl(4, 3, labels = letters[1:4])
  )

  expdat <- data.frame(
    event.a = rep(c(FALSE, TRUE, FALSE), c(0, 3, 9)),
    event.b = rep(c(FALSE, TRUE, FALSE), c(3, 3, 6)),
    event.c = rep(c(FALSE, TRUE, FALSE), c(6, 3, 3)),
    event.d = rep(c(FALSE, TRUE, FALSE), c(9, 3, 0)),
    time = 1:12
  )

  expect_equal(rpaf::expand_events(testdat, "event"), expdat)
})
