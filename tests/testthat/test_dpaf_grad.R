library(testthat)

context("Gradients")
library(rpaf)


# Example data (from whiteboard work) -------------------------------------
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
hz_d <- c(1,2,1,2,4,5,4,5)
hz_m <- c(2,3,2,3,6,7,6,7)
sv_d <- exp(-hz_d)
sv_m <- exp(-hz_m)


# Tests -------------------------------------------------------------------

test_that("Morbidity gradients", {
  grad <- dpaf_grad(z, hz_d, hz_m, sv_d, sv_m, breaks, ID, PERIOD)

  # expected values
  eghz_d <- matrix(c(1,0,0, 0,2,0, 1,0,0, 0,2,0, 4,0,4, 0,5,5, 4,0,4, 0,5,5),
                   ncol = 3, byrow = TRUE)
  eghz_m <- matrix(c(2,0,0, 0,3,0, 2,0,0, 0,3,0, 6,0,6, 0,7,7, 6,0,6, 0,7,7),
                   ncol = 3, byrow = TRUE)
  ecs_ghz_d <- matrix(c(1,0,0, 1,2,0, 1,0,0, 1,2,0,
                        4,0,4, 4,5,9, 4,0,4, 4,5,9),
                      ncol = 3, byrow = TRUE)
  ecs_ghz_m <- matrix(c(2,0,0, 2,3,0, 2,0,0, 2,3,0,
                        6,0,6, 6,7,13, 6,0,6, 6,7,13),
                      ncol = 3, byrow = TRUE)
  egsv_d <- -exp(-hz_d) * ecs_ghz_d
  egsv_m <- -exp(-hz_m) * ecs_ghz_m

  expect_equal(grad$hazard.disease, eghz_d)
  expect_equal(grad$hazard.mortality, eghz_m)
  expect_equal(grad$survival.disease, egsv_d)
  expect_equal(grad$survival.mortality, egsv_m)

  # remains to test final grad-morbidity
})