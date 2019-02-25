library(testthat)

context("Disease PAF Gradients")
library(rpaf)


# Example data (from whiteboard work) -------------------------------------

# Tests -------------------------------------------------------------------

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
  hz_d <- c(1,2,1,2,4,5,4,5)
  hz_m <- c(2,3,2,3,6,7,6,7)
  lambda <- list("primary" = hz_d, "secondary" = hz_m)

  sv_d <- exp(-hz_d)
  sv_m <- exp(-hz_m)
  S <- list("primary" = sv_d, "secondary" = sv_m)

  # grad <- dpaf_grad(z, hz_d, hz_m, sv_d, sv_m, breaks, ID, PERIOD)

  # expected values
  e_glambda_d <- matrix(c(1,0,0, 0,2,0, 1,0,0, 0,2,0, 4,0,4, 0,5,5, 4,0,4, 0,5,5),
                        ncol = 3, byrow = TRUE)
  e_glambda_m <- matrix(c(2,0,0, 0,3,0, 2,0,0, 0,3,0, 6,0,6, 0,7,7, 6,0,6, 0,7,7),
                   ncol = 3, byrow = TRUE)
  ecs_ghz_d <- matrix(c(1,0,0, 1,2,0, 1,0,0, 1,2,0,
                        4,0,4, 4,5,9, 4,0,4, 4,5,9),
                      ncol = 3, byrow = TRUE)
  ecs_ghz_m <- matrix(c(2,0,0, 2,3,0, 2,0,0, 2,3,0,
                        6,0,6, 6,7,13, 6,0,6, 6,7,13),
                      ncol = 3, byrow = TRUE)
  e_gS_d <- -exp(-hz_d) * ecs_ghz_d
  e_gS_m <- -exp(-hz_m) * ecs_ghz_m

  glambda <- rpaf:::dpaf_grad_lambda(z, lambda)
  expect_equal(glambda, list("primary" = e_glambda_d, "secondary" = e_glambda_m))

  gS <- rpaf:::dpaf_grad_S(glambda, S, ID, PERIOD, diff(breaks))
  expect_equal(gS, list("primary" = e_gS_d, "secondary" = e_gS_m))
})

test_that("morbidity gradients", {
  z <- diag(3) # 3x3 identity matrix
  id <- rep(1, 3)
  period <- gl(3, 1)
  dt <- c(1,1,1)

  S <- list("primary" = c(0.95, 0.85, 0.70), "secondary" = c(0.99, 0.95, 0.90))
  lambda <- list("primary" = log(c(1/0.95, 0.95/0.85, 0.85/0.70)),
                 "secondary" = log(c(1/0.99, 0.99/0.95, 0.95/0.90)))

  stopifnot(all.equal(S, rpaf:::dpaf_S(lambda, id, period, dt)))

  glambda <- list("primary" = diag(lambda$primary),
                  "secondary" = diag(lambda$secondary))

  stopifnot(all.equal(glambda, rpaf:::dpaf_grad_lambda(z, lambda)))

  gS <- list("primary" = -S$primary %o% lambda$primary,
             "secondary" = -S$secondary %o% lambda$secondary)
  gS <- lapply(gS, function(mat) {mat[upper.tri(mat)] <- 0; mat})

  stopifnot(all.equal(gS, rpaf:::dpaf_grad_S(glambda, S, id, period, dt)))

  # product of survivals
  Sp <- c(0.9405, 0.8075, 0.6300)

  Sp_gS_d <- (Sp * S[["primary"]]) %o% (-lambda[["primary"]])
  Sp_gS_d[upper.tri(Sp_gS_d)] <- 0
  Sp_gS_m <- (Sp * S[["secondary"]]) %o% (-lambda[["secondary"]])
  Sp_gS_m[upper.tri(Sp_gS_m)] <- 0

  dSp <- -diff(c(1, Sp))
  dSp_gS_d <- apply(Sp_gS_d, 2, function(col) -diff(c(0, col)))
  dSp_gS_m <- apply(Sp_gS_m, 2, function(col) -diff(c(0, col)))
  grad_dSp <- list("primary" = dSp_gS_d, "secondary" = dSp_gS_m)

  dis_prob <- lambda[["primary"]] / do.call(`+`, lambda)

  gdis_prob_d <- lambda[["secondary"]] / do.call(`+`, lambda)**2 * glambda$primary
  gdis_prob_m <- -lambda[["primary"]] / do.call(`+`, lambda)**2 * glambda$secondary

  smnd <- list(
    "primary" = gdis_prob_d * dSp + dis_prob * dSp_gS_d,
    "secondary" = gdis_prob_m * dSp + dis_prob * dSp_gS_m
  )
  rpaf:::dpaf_grad_I(glambda, grad_dSp, lambda, dSp, period)
  expect_equal(rpaf:::dpaf_grad_I(glambda, grad_dSp, lambda, dSp, period),
               smnd)
})
