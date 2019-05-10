context("test-pd_lm")

test_that("invprobit works", {
  xg <- seq(-4, 5)
  expect_equal(invprobit(xg, rho=-1, zeta=2), pnorm(xg, -1, 2, lower.tail = TRUE))
  expect_equal(invprobit(xg, rho=-1, zeta=-2), pnorm(xg, -1, 2, lower.tail = FALSE))
})


test_that("all observed gives equal results", {
  X <- cbind(1,
             rep(c(0.5, -0.5), each=3))
  y <- c(17.6, 17.5, 18.4, 15, 17.5, 16)

  fit <- pd_lm(y ~ X - 1, dropout_curve_position = NA, dropout_curve_scale = NA)
  expect_equal(fit$coefficients, coef(lm(y ~ X - 1)))
  expect_equal(unname(fit$coef_variance_matrix), summary(lm(y ~ X - 1))$sigma^2 * (length(y) - ncol(X)) / length(y) *
                 solve(t(X) %*% X))
})

test_that("some missing gives equal results", {

  X <- cbind(1,
             rep(c(0.5, -0.5), each=3))
  y <- c(17.6, NA, NA, 15, 17.5, 16)

  fit1 <- pd_lm(y ~ X - 1, dropout_curve_position = 18, dropout_curve_scale = -1, method = "analytic_hessian")
  fit2 <- pd_lm(y ~ X - 1, dropout_curve_position = 18, dropout_curve_scale = -1, method = "analytic_grad")
  fit3 <- pd_lm(y ~ X - 1, dropout_curve_position = 18, dropout_curve_scale = -1, method = "numeric")
  expect_equal(fit1, fit2, tolerance = 1e-2)
  expect_equal(fit3, fit2, tolerance = 1e-2)
})


test_that("pd_lm works with moderaton", {

  X <- matrix(rep(c(0.5, -0.5), each=3), ncol=1)
  y <- c(17.6, 17.5, 18.4, NA, 17.5, NA)
  t_rho <- rep(18, 6)
  t_zeta <- rep(-0.8, 6)
  t_mu0 <- 20
  t_sigma20 <- 5
  t_tau20 <- 0.04
  t_df0 <- 1.3

  fit1 <- pd_lm(y ~ X, dropout_curve_position = t_rho, dropout_curve_scale = t_zeta,
                location_prior_mean = t_mu0, location_prior_scale = t_sigma20,
                variance_prior_scale = t_tau20, variance_prior_df = t_df0,
                method = "analytic_hessian")
  fit2 <- pd_lm(y ~ X, dropout_curve_position = t_rho, dropout_curve_scale = t_zeta,
                location_prior_mean = t_mu0, location_prior_scale = t_sigma20,
                variance_prior_scale = t_tau20, variance_prior_df = t_df0,
                method = "analytic_grad")
  fit3 <- pd_lm(y ~ X, dropout_curve_position = t_rho, dropout_curve_scale = t_zeta,
                location_prior_mean = t_mu0, location_prior_scale = t_sigma20,
                variance_prior_scale = t_tau20, variance_prior_df = t_df0,
                method = "numeric")
  expect_equal(fit1, fit2, tolerance = 1e-2)
  expect_equal(fit3, fit2, tolerance = 1e-1)

})



test_that("pd_lm has the right names", {
  col_data <- data.frame(f1 = factor(rep(LETTERS[1:5],each=2), levels = LETTERS[1:7]),
                         f2 = factor(c("Good", "Neutral", "Neutral", "Bad", "Bad",
                                       "Bad", "Good", "Bad", "Neutral", "Bad"),
                                     levels = c("Bad", "Neutral", "Good"), ordered=TRUE),
                         f3 = factor(rep(c("hello", "world"),times=5)),
                         f4 = factor(rep(c("hello", "world", "foo", "bar", "foobar"),times=2)),
                         c1 = sample(rep(c("ABC", "xyz"), each=5)),
                         num = rnorm(10),
                         num2 = rnorm(10),
                         stringsAsFactors = FALSE)

  y <- rnorm(10, mean=10)
  col_data$resp <- y
  lm_coef <- coef(lm(y ~ f4, data=col_data))
  fit <- pd_lm(y ~ f4, data = col_data, dropout_curve_position = NA, dropout_curve_scale = NA)
  fit1.5 <- pd_lm(resp ~ f4, data = col_data, dropout_curve_position = NA, dropout_curve_scale = NA)
  fit2 <- pd_lm(y ~ f4 - 1, data = col_data, dropout_curve_position = NA, dropout_curve_scale = NA)

  expect_equal(unname(lm_coef), unname(fit$coefficients))
  expect_equal(names(fit$coefficients[-1]),
               names(fit2$coefficients[-1]))
  expect_equal(fit$coefficients["Intercept"] + fit$coefficients[-1],
               fit2$coefficients[-1])

})





test_that("derivatives are correct", {

  dropout_curve_position <- rep(8, times=5)
  dropout_curve_scale <- rep(-1, times=5)
  mu0 <- 14
  sigma20 <- 0.3
  df0 <- 2
  tau20 <- 0.05
  location_prior_df <- 4
  moderate_location <- FALSE
  moderate_variance <- FALSE

  X <- cbind(c(1,1,0,0,0),c(0,0,1,1,1))
  y <- c(13, 12, NA, 11, 10)

  Xo <- X[!is.na(y), , drop=FALSE]
  Xm <- X[is.na(y), , drop=FALSE]
  yo <- y[!is.na(y)]
  p <- ncol(X)
  rho <- dropout_curve_position[is.na(y)]
  zeta <- dropout_curve_scale[is.na(y)]

  par <- c(4, 4, 1)
  num_grad <- numDeriv::grad(func = function(par){
    beta <- par[seq_len(ncol(X))]
    sigma2 <- par[p+1]
    if(sigma2 <= 0) return(10000)
    zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
    objective_fnc(y, yo, X, Xm, Xo,
                  beta, sigma2, rho, zetastar, mu0, sigma20, df0, tau20,
                  location_prior_df, moderate_location, moderate_variance)
  }, x = par)

  num_hessian <- numDeriv::hessian(func = function(par){
    beta <- par[seq_len(ncol(X))]
    sigma2 <- par[p+1]
    if(sigma2 <= 0) return(10000)
    zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
    objective_fnc(y, yo, X, Xm, Xo,
                  beta, sigma2, rho, zetastar, mu0, sigma20, df0, tau20,
                  location_prior_df, moderate_location, moderate_variance)
  }, x = par)

  beta <- par[seq_len(ncol(X))]
  sigma2 <- par[p+1]
  zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
  grad_result <- grad_fnc(y, yo, X, Xm, Xo,
                beta, sigma2, rho, zetastar, mu0, sigma20, df0, tau20,
                location_prior_df, moderate_location, moderate_variance)
  hessian_result <- hess_fnc(y, yo, X, Xm, Xo,
                 beta, sigma2, rho, zetastar, mu0, sigma20, df0, tau20,
                 location_prior_df, moderate_location, moderate_variance,
                 1:2, 2)
  expect_equal(num_grad, grad_result, tol = 1e-8)
  expect_equal(num_hessian, hessian_result, tol = 1e-8)








  num_grad2 <- numDeriv::grad(func = function(par){
    beta <- par[seq_len(ncol(X))]
    sigma2 <- par[p+1]
    zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
    objective_fnc(y, yo, X, Xm, Xo,
                  beta, sigma2, rho, zetastar, mu0, sigma20, df0, tau20,
                  location_prior_df, moderate_location=TRUE, moderate_variance=TRUE)
  }, x = par)

  num_hessian2 <- numDeriv::hessian(func = function(par){
    beta <- par[seq_len(ncol(X))]
    sigma2 <- par[p+1]
    zetastar <- zeta * sqrt(1 + sigma2/zeta^2)
    objective_fnc(y, yo, X, Xm, Xo,
                  beta, sigma2, rho, zetastar, mu0, sigma20, df0, tau20,
                  location_prior_df, moderate_location=TRUE, moderate_variance=TRUE)
  }, x = par)
  grad_result2 <- grad_fnc(y, yo, X, Xm, Xo,
                          beta, sigma2, rho, zetastar, mu0, sigma20, df0, tau20,
                          location_prior_df, moderate_location=TRUE, moderate_variance=TRUE)
  hessian_result2 <- hess_fnc(y, yo, X, Xm, Xo,
                             beta, sigma2, rho, zetastar, mu0, sigma20, df0, tau20,
                             location_prior_df, moderate_location=TRUE, moderate_variance=TRUE,
                             1:2, 2)
  expect_equal(num_grad2, grad_result2, tol = 1e-8)
  expect_equal(num_hessian2, hessian_result2, tol = 1e-8)



})





test_that("available information is correctly reflected", {

  # var and loc moderation --> everything fine
  y <- c(NA, NA, NA, NA)
  X <- cbind(c(1,1,0,0), c(0,0,1,1))
  fit_all_mis <- pd_lm(y ~ X - 1, dropout_curve_position = 0, dropout_curve_scale = -1,
        location_prior_mean = 0, location_prior_scale = 1,
        variance_prior_scale = 0.1, variance_prior_df = 2)
  expect_true(all(! is.na(coefficients(fit_all_mis))))
  expect_true(! is.na(fit_all_mis$s2))

  # only location moderation --> At least one condition with one obs
  # Sometimes the information is not enough than all is NA...
  y <- c(NA, NA, NA, NA)
  X <- cbind(c(1,1,0,0), c(0,0,1,1))
  fit_all_mis <- pd_lm(y ~ X - 1, dropout_curve_position = 0, dropout_curve_scale = -1,
                       location_prior_mean = 0, location_prior_scale = 1,
                       location_prior_df = 100)
  expect_true(all(is.na(coefficients(fit_all_mis))))
  expect_true(is.na(fit_all_mis$s2))

  y <- c(5, NA, NA, NA, NA)
  X <- cbind(c(1,1,0,0, 0), c(0,0,1,1,0), c(0,0,0,0,1))
  fit_all_but_on_mis <- pd_lm(y ~ X - 1, dropout_curve_position = 0, dropout_curve_scale = -1,
                              location_prior_mean = 0, location_prior_scale = 1)
  expect_true(all(!is.na(coefficients(fit_all_but_on_mis))))
  expect_true(!is.na(fit_all_but_on_mis$s2))

  # var moderation --> Everything is fine
  # but if there is no obs for a cond, than beta = NA
  y <- c(NA, NA, NA, NA)
  X <- cbind(c(1,1,0,0), c(0,0,1,1))
  fit_all_mis <- pd_lm(y ~ X - 1, dropout_curve_position = 0, dropout_curve_scale = -1,
                       variance_prior_scale = 0.1, variance_prior_df = 2)
  expect_true(all(is.na(coefficients(fit_all_mis))))
  expect_true(is.na(fit_all_mis$s2))

  y <- c(5, NA, NA, NA, NA)
  X <- cbind(c(1,1,1,1, 1), c(0,0,1,1,0), c(0,0,0,0,1))
  fit_all_but_on_mis <- pd_lm(y ~ X - 1, dropout_curve_position = 0, dropout_curve_scale = -1,
                              variance_prior_scale = 0.1, variance_prior_df = 2)
  expect_true(!is.na(coefficients(fit_all_but_on_mis)[1]))
  expect_true(all(is.na(coefficients(fit_all_but_on_mis)[2:3])))
  expect_true(! is.na(fit_all_but_on_mis$s2))

  # No moderation. At least one observation.
  # but if there is no obs for a cond, than beta = -Inf
  y <- c(NA, NA, NA, NA)
  X <- cbind(c(1,1,0,0), c(0,0,1,1))
  fit_all_mis <- pd_lm(y ~ X - 1, dropout_curve_position = 0, dropout_curve_scale = -1)
  expect_true(all(is.na(coefficients(fit_all_mis))))
  expect_true(is.na(fit_all_mis$s2))

  y <- c(5, NA, NA, NA, NA)
  X <- cbind(c(1,1,0,0, 0), c(0,0,1,1,0), c(0,0,0,0,1))
  fit_all_but_on_mis <- pd_lm(y ~ X - 1, dropout_curve_position = 0, dropout_curve_scale = -1)
  expect_true(!is.na(coefficients(fit_all_but_on_mis)[1]))
  expect_true(all(is.na(coefficients(fit_all_but_on_mis)[2:3])))
  expect_true(!is.na(fit_all_but_on_mis$s2))

})



test_that("There are more values than parameters", {
  y <- c(12, 13)
  X <- matrix(c(1,0, 0, 1), ncol=2, nrow=2)

  expect_silent(pd_lm(y ~ X - 1, dropout_curve_position = NA, dropout_curve_scale = NA,
                      location_prior_mean = 0, location_prior_scale = 1,
                      variance_prior_scale = 0.1, variance_prior_df = 2))
  expect_silent(pd_lm(y ~ X - 1, dropout_curve_position = NA, dropout_curve_scale = NA,
                      variance_prior_scale = 0.1, variance_prior_df = 2))

  expect_silent(pd_lm(y ~ X - 1, dropout_curve_position = NA, dropout_curve_scale = NA,
                      location_prior_mean = 0, location_prior_scale = 1))

  expect_error(pd_lm(y ~ X - 1, dropout_curve_position = NA, dropout_curve_scale = NA))


})







test_that("parametrization does not influence result", {

  y <- c(23,24.1, 25.5, NA, NA, NA)
  X1 <- matrix(c(1,1,1,1,1,1,   0,0,0,1,1,1), ncol=2)
  X2 <- matrix(c(1,1,1,0,0,0,   0,0,0,1,1,1), ncol=2)

  rho <- 22.2
  zeta <- -0.3

  f1 <- pd_lm(y ~ X1 - 1, dropout_curve_position = rho, dropout_curve_scale = zeta)
  f2 <- pd_lm(y ~ X2 - 1, dropout_curve_position = rho, dropout_curve_scale = zeta)
  f1$coefficients
  f2$coefficients


  # This one fails because of problematic initialization
  expect_equal(unname(f1$coefficients),
               unname(c(f2$coefficients[1], f2$coefficients[2] - f2$coefficients[1])),
               tolerance = 1e-2)
  expect_equal(f1$s2, f2$s2, tolerance=1e-6)
})






