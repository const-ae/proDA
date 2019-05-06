context("test-calculate_distance")


test_that("distance_sq method works", {
  set.seed(1)
  means_x <- rnorm(50, mean=20, sd=2)
  means_y <- rnorm(50, mean=17, sd=2)
  var_x <- rep(3.1^2, length(means_x))
  var_y <- rep(0.3^2, length(means_x))

  X <- sapply(seq_along(means_x), function(idx){
    rnorm(500, mean=means_x[idx], sd=sqrt(var_x[idx]))
  })
  Y <- sapply(seq_along(means_y), function(idx){
    rnorm(500, mean=means_y[idx], sd=sqrt(var_y[idx]))
  })
  res <- distance_sq(mu1 = means_x, sigma1 = var_x,
                     mu2 = means_y, sigma2 = var_y)

  dists_sq <- rowSums((X - Y)^2)
  dists <- sqrt(dists_sq)
  expect_equal(mean(dists_sq) / res$mean, 1, tolerance = 1e-1)
  expect_equal(sd(dists_sq) / sqrt(res$var), 1, tolerance = 1e-1)
  expect_equal(mean(dists) / sqrt(res$mean), 1, tolerance = 1e-1)
  expect_equal(sd(dists) / sqrt(res$var * 1/(4 * res$mean)), 1, tolerance = 1e-1)
  # xg <- seq(min(dists_sq), max(dists_sq), length.out= 1001)
  # hist(dists_sq, prob= TRUE)
  # lines(xg, dnorm(xg, mean(dists_sq), sd(dists_sq)), col="black" )
  # lines(xg, dnorm(xg, res$mean, sqrt(res$var)), col="red" )
  # xg2 <- seq(min(sqrt(dists_sq)), max(sqrt(dists_sq)), length.out= 1001)
  # hist(sqrt(dists_sq), prob= TRUE)
  # lines(xg2, dnorm(xg2, mean(sqrt(dists_sq)), sd(sqrt(dists_sq))), col="black" )
  # lines(xg2, dnorm(xg2, sqrt(res$mean), sqrt(res$var * 1/(4 * res$mean)) ), col="red" )
})



test_that("dist_approx works", {

  se <- generate_synthetic_data(n_proteins = 200, n_conditions = 4,
                                dropout_curve_position = 17.2,
                                location_prior_scale = 1,
                                n_replicates = 3,
                                return_summarized_experiment = TRUE)
  fit <- proDA(se, c(rep(LETTERS[1:4], each=3)),
               verbose=TRUE)

  da <- dist_approx(fit)
  da2 <- dist_approx(fit, by_sample = FALSE)

  expect_gt(cor(c(dist(t(assay(se, 2)))), c(da$mean)), 0.7)
  expect_gt(cor(c(dist(assay(se, 2))), c(da2$mean)), 0.7)

})

