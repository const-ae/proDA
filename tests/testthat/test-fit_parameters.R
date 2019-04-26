context("test-fit_parameters")

test_that("fit_parameters works", {
  set.seed(1)
  data <- matrix(rnorm(1000 * 5), nrow=1000, ncol=5)

  fit <- fit_parameters(data, location_prior_df = 100, verbose=TRUE)
  fit
  test_res <- test_diff(fit, "Intercept")
  test_res

  var(rowMeans(fit$coefficients %*% t(fit$design)))
  fit$coefficients %*% t(fit$design)
  plot(rowMeans(data), test_res$avg_abundance); abline(0,1)
})
