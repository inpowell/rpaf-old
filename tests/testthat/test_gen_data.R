library(testthat)

context("test gen_data function")

library(rpaf)

df <- data.frame(
  id = 1:6,
  tmt = gl(2, 3, labels = c("Tmt", "Ctl")),
  ind = rep(c(TRUE, FALSE), 3),
  time = 5 + (1:6) %% 3
)

test_that("unordered breaks are corrected with warning", {
  expect_warning(paf_data <- gen_data(df, "id", ft_breaks = c(0,2,6,4),
                                      death_ind = "ind", death_time = "time",
                                      variables = "tmt"),
                 "Time breaks are not in order; reordering them.")
  expect_equal(paf_data$breaks, c(0,2,4,6))
})
