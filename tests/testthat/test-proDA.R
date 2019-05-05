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
  pd_res <- test_diff(pd_fit, "A")

  lim_fit <- limma::lmFit(data, pd_fit$design)
  lim_fit <- limma::eBayes(lim_fit)
  lim_res <- limma::topTable(lim_fit, "A", sort.by = "none", number = 100)

  expect_equal(unname(lim_fit$coefficients),
               unname(pd_fit$coefficients), tolerance = 1e-6)
  expect_equal(lim_fit$s2.post, pd_fit$feature_parameters$s2, tolerance = 1e-3)

  expect_equal(lim_res$logFC, pd_res$diff, tolerance = 1e-6)
  expect_equal(lim_res$AveExpr, pd_res$avg_abundance, tolerance = 1e-6)
  expect_equal(lim_res$t, pd_res$t_statistic, tolerance = 1e-2)
  expect_equal(lim_res$P.Value, pd_res$pval, tolerance = 1e-2)
  expect_equal(lim_res$adj.P.Val, pd_res$adj_pval, tolerance = 1e-3)
})




test_that("proDA works with missing values", {
  skip("Very slow")
  # This is only here to fail if there are errors in the code
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

  expect_true(pd_fit1$convergence$successful)
  expect_true(pd_fit2$convergence$successful)
  expect_true(pd_fit3$convergence$successful)
  expect_true(pd_fit4$convergence$successful)
})

test_that("proDA works with many missing values", {
  skip("Very slow")
  # This is only here to fail if there are errors in the code
  set.seed(1)
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


  expect_true(pd_fit1$convergence$successful)
  expect_true(pd_fit2$convergence$successful)
  expect_true(pd_fit3$convergence$successful)
  expect_true(pd_fit4$convergence$successful)
})



test_that("predict works", {
  set.seed(1)
  data <- matrix(rnorm(100 * 5), nrow=100, ncol=5)
  colnames(data) <- paste0("sample_", 1:5)
  rownames(data) <- paste0("protein_", 1:100)

  fit <- proDA(data, moderate_location = FALSE, verbose=TRUE)

  Pred <- predict(fit, type="response")
  pred_fit <- predict(fit, type="feature_parameters")

  expect_equal(fit$hyper_parameters, pred_fit$hyper_parameters)
  expect_equal(fit$colData, pred_fit$colData)

  # Because there was no 'newdata'
  expect_equal(fit, pred_fit)

  newdata <-  matrix(rnorm(100 * 5), nrow=100, ncol=5)
  colnames(newdata) <- paste0("sample_", 1:5)
  rownames(newdata) <- paste0("protein_", 1:100)

  newfit <- predict(fit, newdata, type="feature_parameters")
  expect_equal(fit$hyper_parameters, newfit$hyper_parameters)
  expect_equal(fit$colData, newfit$colData)


  combfit <- rbind(fit, newfit)
  expect_equal(fit$abundances, combfit[1:100, ]$abundances)
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
  expect_true(fit$convergence$successful)


})





