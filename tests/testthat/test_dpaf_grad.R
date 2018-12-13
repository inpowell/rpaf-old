library(testthat)

context("Gradients")
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
  hz <- cbind(disease = hz_d, mortality = hz_m)

  sv_d <- exp(-hz_d)
  sv_m <- exp(-hz_m)
  sv <- cbind(disease = sv_d, mortality = sv_m)

  # grad <- dpaf_grad(z, hz_d, hz_m, sv_d, sv_m, breaks, ID, PERIOD)

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

  ghz <- rpaf:::dpaf_ghz(z, hz)
  expect_equal(ghz, list(disease = eghz_d, mortality = eghz_m))

  gsv <- rpaf:::dpaf_gsv(ghz, sv, ID, PERIOD, diff(breaks))
  expect_equal(gsv, list(disease = egsv_d, mortality = egsv_m))
})

test_that("morbidity gradients", {
  z <- diag(3) # 3x3 identity matrix
  id <- rep(1, 3)
  period <- gl(3, 1)
  dt <- c(1,1,1)

  sv <- cbind(disease = c(0.95, 0.85, 0.70), mortality = c(0.99, 0.95, 0.90))
  hz <- cbind(disease = log(c(1/0.95, 0.95/0.85, 0.85/0.70)),
              mortality = log(c(1/0.99, 0.99/0.95, 0.95/0.90)))

  stopifnot(all.equal(sv, rpaf:::dpaf_sv(hz, id, period, dt)))

  ghz <- list(disease = diag(hz[,1]), mortality = diag(hz[,2]))

  stopifnot(all.equal(ghz, rpaf:::dpaf_ghz(z, hz)))

  gsv <- list(disease = -sv[,1] %o% hz[,1], mortality = -sv[,2] %o% hz[,2])
  gsv <- lapply(gsv, function(mat) {mat[upper.tri(mat)] <- 0; mat})

  stopifnot(all.equal(gsv, rpaf:::dpaf_gsv(ghz, sv, id, period, dt)))

  # product of survivals
  svp <- c(0.9405, 0.8075, 0.6300)

  svp_gsv_d <- (svp * sv[,"disease"]) %o% (-hz[,"disease"])
  svp_gsv_d[upper.tri(svp_gsv_d)] <- 0
  svp_gsv_m <- (svp * sv[,"mortality"]) %o% (-hz[,"mortality"])
  svp_gsv_m[upper.tri(svp_gsv_m)] <- 0

  dsvp <- -diff(c(1, svp))
  dsvp_gsv_d <- apply(svp_gsv_d, 2, function(col) -diff(c(0, col)))
  dsvp_gsv_m <- apply(svp_gsv_m, 2, function(col) -diff(c(0, col)))

  dis_prob <- hz[,"disease"] / rowSums(hz)

  gdis_prob_d <- hz[,"mortality"] / rowSums(hz)**2 * ghz$disease
  gdis_prob_m <- -hz[,"disease"] / rowSums(hz)**2 * ghz$mortality

  smnd <- list(
    disease = gdis_prob_d * dsvp + dis_prob * dsvp_gsv_d,
    mortality = gdis_prob_m * dsvp + dis_prob * dsvp_gsv_m
  )

  expect_equal(rpaf:::dpaf_gi(ghz, gsv, hz, sv, id, period),
               lapply(smnd, colSums))
})
