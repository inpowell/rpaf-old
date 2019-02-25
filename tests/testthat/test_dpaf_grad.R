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
  hz_1 <- c(1,2,1,2,4,5,4,5)
  hz_2 <- c(2,3,2,3,6,7,6,7)
  lambda <- list("primary" = hz_1, "secondary" = hz_2)

  sv_1 <- exp(-hz_1)
  sv_2 <- exp(-hz_2)
  S <- list("primary" = sv_1, "secondary" = sv_2)

  # grad <- dpaf_grad(z, hz_1, hz_2, sv_1, sv_2, breaks, ID, PERIOD)

  # expected values
  e_glambda_1 <- matrix(c(1,0,0, 0,2,0, 1,0,0, 0,2,0, 4,0,4, 0,5,5, 4,0,4, 0,5,5),
                        ncol = 3, byrow = TRUE)
  e_glambda_2 <- matrix(c(2,0,0, 0,3,0, 2,0,0, 0,3,0, 6,0,6, 0,7,7, 6,0,6, 0,7,7),
                   ncol = 3, byrow = TRUE)
  ecs_ghz_1 <- matrix(c(1,0,0, 1,2,0, 1,0,0, 1,2,0,
                        4,0,4, 4,5,9, 4,0,4, 4,5,9),
                      ncol = 3, byrow = TRUE)
  ecs_ghz_2 <- matrix(c(2,0,0, 2,3,0, 2,0,0, 2,3,0,
                        6,0,6, 6,7,13, 6,0,6, 6,7,13),
                      ncol = 3, byrow = TRUE)
  e_gS_1 <- -exp(-hz_1) * ecs_ghz_1
  e_gS_2 <- -exp(-hz_2) * ecs_ghz_2

  glambda <- rpaf:::dpaf_grad_lambda(z, lambda)
  expect_equal(glambda, list("primary" = e_glambda_1, "secondary" = e_glambda_2))

  gS <- rpaf:::dpaf_grad_S(glambda, S, ID, PERIOD, diff(breaks))
  expect_equal(gS, list("primary" = e_gS_1, "secondary" = e_gS_2))
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

  Sp_gS_1 <- (Sp * S[["primary"]]) %o% (-lambda[["primary"]])
  Sp_gS_1[upper.tri(Sp_gS_1)] <- 0
  Sp_gS_2 <- (Sp * S[["secondary"]]) %o% (-lambda[["secondary"]])
  Sp_gS_2[upper.tri(Sp_gS_2)] <- 0

  dSp <- -diff(c(1, Sp))
  dSp_gS_1 <- apply(Sp_gS_1, 2, function(col) -diff(c(0, col)))
  dSp_gS_2 <- apply(Sp_gS_2, 2, function(col) -diff(c(0, col)))
  grad_dSp <- list("primary" = dSp_gS_1, "secondary" = dSp_gS_2)

  dis_prob <- lambda[["primary"]] / do.call(`+`, lambda)

  gdis_prob_1 <- lambda[["secondary"]] / do.call(`+`, lambda)**2 * glambda$primary
  gdis_prob_2 <- -lambda[["primary"]] / do.call(`+`, lambda)**2 * glambda$secondary

  smnd <- list(
    "primary" = gdis_prob_1 * dSp + dis_prob * dSp_gS_1,
    "secondary" = gdis_prob_2 * dSp + dis_prob * dSp_gS_2
  )
  rpaf:::dpaf_grad_I(glambda, grad_dSp, lambda, dSp, period)
  expect_equal(rpaf:::dpaf_grad_I(glambda, grad_dSp, lambda, dSp, period),
               smnd)
})
