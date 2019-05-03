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

  pd_fit1 <- proDA(data, c("A", "A", "B", "B", "B"),
                   moderate_location = TRUE, verbose=TRUE)
  pd_fit2 <- proDA(data, c("A", "A", "B", "B", "B"),
                   moderate_location = FALSE, verbose=TRUE)
  pd_fit3 <- proDA(data, c("A", "A", "B", "B", "B"),
                   moderate_variance = FALSE, verbose=TRUE)
  pd_fit4 <- proDA(data, c("A", "A", "B", "B", "B"),
                   moderate_variance = FALSE, moderate_location = FALSE,
                   verbose=TRUE)
})

test_that("proDA works with many missing values", {
  set.seed(2)
  data <- matrix(rnorm(100 * 5), nrow=100, ncol=5)
  data[sample(seq_len(100 * 5), 350)] <- NA

  pd_fit1 <- proDA(data, c("A", "A", "B", "B", "B"),
                   moderate_location = TRUE, verbose=TRUE)
  pd_fit2 <- proDA(data, c("A", "A", "B", "B", "B"),
                  moderate_location = FALSE, verbose=TRUE)
  pd_fit3 <- proDA(data, c("A", "A", "B", "B", "B"),
                  moderate_variance = FALSE, verbose=TRUE)
  pd_fit4 <- proDA(data, c("A", "A", "B", "B", "B"),
                   moderate_variance = FALSE, moderate_location = FALSE,
                   verbose=TRUE)


  head(test_diff(pd_fit1, A - B), n=5)
  head(test_diff(pd_fit4, A - B), n=5)
})



test_that("predict works", {
  set.seed(1)
  data <- matrix(rnorm(100 * 5), nrow=100, ncol=5)
  colnames(data) <- paste0("sample_", 1:5)
  rownames(data) <- paste0("protein_", 1:100)

  fit <- proDA(data, moderate_location = FALSE, verbose=TRUE)
  test_res <- test_diff(fit, "Intercept")

  expect_equal(rowMeans(data), test_res$avg_abundance)

  Pred <- predict(fit, type="response")
  predict(fit, type="feature_parameters")
  newdata <-  matrix(rnorm(100 * 5), nrow=100, ncol=5)
  colnames(newdata) <- paste0("sample_", 1:5)
  rownames(newdata) <- paste0("protein_", 1:100)

  newfit <- predict(fit, newdata, type="feature_parameters")
  combfit <- rbind(fit, newfit)

  combfit[c(1,101), ]$abundances
  predict(fit, newdata)
})



test_that("fit works with SummarizedExperiment object", {

  data <- matrix(rnorm(22 * 50, mean=20), ncol=22, nrow=50)
  data[invprobit(data, 18, -1) > runif(22 * 50)] <- NA
  colnames(data) <- paste0("sample_", LETTERS[seq_len(ncol(data))])
  rownames(data) <- paste0("protein", 1:50)
  annot_df <- data.frame(Name = paste0("sample_", LETTERS[seq_len(ncol(data))]),
                         cond = c(rep(c("A", "B", "C"), times=7), "D"),
                         num = runif(ncol(data)))
  se <- SummarizedExperiment(data, colData = annot_df)

  fit <- proDA(se, ~ cond + num, verbose=TRUE)
  fit$design
  test_diff(fit, num)


})


test_that("t-test works", {
  set.seed(1)
  data <- matrix(rnorm(100 * 5), nrow=100, ncol=5)
  data[sample(seq_len(100 * 5), 150)] <- NA
  colnames(data) <- paste0("sample_", 1:5)
  rownames(data) <- paste0("protein_", 1:100)


  pd_row_t_test(X = data[,1:2], Y=data[,3:5])

})





test_that("parameter recovery works", {

  data <- generate_synthetic_data(n_proteins = 500,
                                  dropout_curve_scale = -0.4)
  hist(data$Z, col=ggplot2::alpha("green", 0.3))
  hist(data$Y, add=TRUE, col=ggplot2::alpha("red", 0.2))
  fit <- proDA(data$Y, design = data$groups, verbose=TRUE)


  # dat <- generate_synthetic_data(n_proteins = 100, frac_changed = 0)$Z
  dat <- matrix(rnorm(100 * 6, mean=20, sd=3), nrow=100)
  dat[sample(seq_len(prod(dim(dat))), prod(dim(dat)) * 0.6)] <- NA
  fit2 <- proDA(dat, data$groups,
                moderate_location = FALSE, moderate_variance = FALSE,
                verbose=TRUE)
  fit2

})


