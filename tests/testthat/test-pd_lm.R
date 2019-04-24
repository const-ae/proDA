context("test-pd_lm")

test_that("invprobit works", {
  xg <- seq(-4, 5)
  expect_equal(invprobit(xg, rho=-1, zeta=2), pnorm(xg, -1, 2, lower.tail = TRUE))
  expect_equal(invprobit(xg, rho=-1, zeta=-2), pnorm(xg, -1, 2, lower.tail = FALSE))
})



test_that("pd_lm works with missing values", {

  X <- cbind(1,
             rep(c(0.5, -0.5), each=3))
  y <- c(17.6, 17.5, 18.4, NA, 17.5, NA)
  t_rho <- rep(18, 6)
  t_zeta <- rep(-0.8, 6)
  t_mu0 <- 20
  t_sigma20 <- 5
  t_tau20 <- 0.04
  t_df0 <- 1.3
  fit <- pd_lm(y, X,
     rho = t_rho, zeta = t_zeta, mu0 = t_mu0, sigma20 = t_sigma20,
     tau20 = t_tau20, df0 = t_df0, method="analytic")

  lm(y ~ X)
  fit

})


test_that("pd_lm_unreg works without missing values", {

  X <- cbind(1,
             rep(c(0.5, -0.5), each=3))
  y <- c(17.6, 17.5, 18.4, 17.2, 17.5, 17.8)
  t_rho <- rep(18, 6)
  t_zeta <- rep(-0.8, 6)

  fit <- pd_lm_unreg(y, X,
    rho = t_rho, zeta = t_zeta,  method="analytic")
  lin_m <- lm(y ~ X)
  expect_equal(length(y), fit$n_approx)
  expect_equal(summary(lin_m)$sigma, sqrt(fit$s2))
  expect_equal(sum((y - predict(lin_m))^2), fit$rss)

})







