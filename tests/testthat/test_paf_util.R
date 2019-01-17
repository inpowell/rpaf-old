library(testthat)

context("PAF generic utitlities")

library(rpaf)

test_that("Period sorting", {
  id <- gl(3, 3)
  period_sorted <- gl(3, 1, 9)
  period_unsorted <- period_sorted[c(1:7, 9, 8)]

  expect_false(rpaf:::is_period_unsorted(id, period_sorted))
  expect_true(rpaf:::is_period_unsorted(id, period_unsorted))
})
