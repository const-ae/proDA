context("test-proDA")

test_that("proDA works", {
  set.seed(1)
  data <- matrix(rnorm(100 * 5), nrow=100, ncol=5)

  fit <- proDA(data, moderate_location = FALSE, verbose=TRUE)
  test_res <- test_diff(fit, "Intercept")

  expect_equal(rowMeans(data), test_res$avg_abundance)


})


test_that("proDA works as good as limma with observed values", {
  set.seed(1)
  data <- matrix(rnorm(100 * 5), nrow=100, ncol=5)

  pd_fit <- proDA(data, c("A", "A", "B", "B", "B"),
               moderate_location = FALSE, verbose=TRUE)
  pd_res <- test_diff(fit, "A")

  lim_fit <- limma::lmFit(data, pd_fit$design)
  lim_fit <- limma::eBayes(lim_fit)
  lim_res <- limma::topTable(lim_fit, "A", sort.by = "none", number = 100)

  expect_equal(unname(lim_fit$coefficients),
               unname(pd_fit$coefficients))
  expect_equal(lim_fit$s2.post, pd_fit$feature_parameters$s2, tolerance = 1e-3)

  expect_equal(lim_res$logFC, pd_res$diff)
  expect_equal(lim_res$AveExpr, pd_res$avg_abundance)
  expect_equal(lim_res$t, pd_res$t_statistic, tolerance = 1e-3)
  expect_equal(lim_res$P.Value, pd_res$pval, tolerance = 1e-3)
  expect_equal(lim_res$adj.P.Val, pd_res$adj_pval, tolerance = 1e-3)
})


test_that("dropout curve idenfication works", {
  n_prots <- 1000
  t_Y <- matrix(rnorm(n_prots), ncol=1)
  Y <- t_Y
  Y[sample(seq_len(n_prots), 75)] <- NA
  X <- matrix(1, ncol=1, nrow=1)
  Pred <- matrix(rnorm(n_prots, mean=t_Y, sd=0.1), ncol=1)
  s2 <- rep(0.1, n_prots)

  dropout_curves(Y, X, Pred, s2)

})


test_that("proDA works with missing values", {
  set.seed(1)
  data <- matrix(rnorm(100 * 5), nrow=100, ncol=5)
  data[sample(seq_len(100 * 5), 150)] <- NA

  pd_fit <- proDA(data, c("A", "A", "B", "B", "B"),
                  moderate_location = FALSE, verbose=TRUE)
})

test_that("proDA works with many missing values", {
  set.seed(1)
  data <- matrix(rnorm(100 * 5), nrow=100, ncol=5)
  data[sample(seq_len(100 * 5), 150)] <- NA

  pd_fit <- proDA(data, c("A", "A", "B", "B", "B"),
                  moderate_location = FALSE, verbose=TRUE)
  pd_fit2 <- proDA(data, c("A", "A", "B", "B", "B"),
                  moderate_location = TRUE, verbose=TRUE)
})






