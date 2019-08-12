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
  # data <- matrix(rnorm(100 * 5), nrow=100, ncol=5)
  data <- generate_synthetic_data(100, n_replicates = c(2,3), frac_changed = 0)$Z

  pd_fit <- proDA(data, c("A", "A", "B", "B", "B"),
               moderate_location = FALSE, verbose=TRUE)
  pd_res <- test_diff(pd_fit, A - B)

  lim_fit <- limma::lmFit(data, pd_fit$design)
  lim_fit <- limma::contrasts.fit(lim_fit, limma::makeContrasts(A - B, levels = result_names(pd_fit)))
  lim_fit <- limma::eBayes(lim_fit)
  lim_res <- limma::topTable(lim_fit, sort.by = "none", number = 100)

  expect_equal(unname(c(lim_fit$coefficients)),
               unname(c(pd_fit$coefficients[,1] - pd_fit$coefficients[,2])), tolerance = 1e-6)
  expect_gt(cor(unname(log(lim_fit$s2.post)), unname(log(pd_fit$feature_parameters$s2))), 0.99)

  expect_equal(lim_res$logFC, pd_res$diff, tolerance = 1e-6)
  expect_equal(lim_res$AveExpr, pd_res$avg_abundance, tolerance = 1e-6)
  expect_gt(cor(lim_res$t, pd_res$t_statistic), 0.99)
  expect_gt(cor(lim_res$P.Value, pd_res$pval), 0.99)
})




test_that("proDA works with missing values", {
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

# test_that("fit works with MSnSet", {
#   library(MSnbase)
#
#   syn_data <- generate_synthetic_data(n_proteins = 20)
#
#   annotation_df <- data.frame(group = syn_data$groups)
#   rownames(annotation_df) <- colnames(syn_data$Y)
#
#   fData <- data.frame(true_mean1 = syn_data$t_mu[,1],
#                       true_mean2 = syn_data$t_mu[,2],
#                       changed= syn_data$changed)
#   rownames(fData) <- rownames(syn_data$Y)
#   ms <- MSnSet(syn_data$Y, pData=AnnotatedDataFrame(annotation_df),
#                fData=AnnotatedDataFrame(fData))
#
#   expect_equal(median_normalization(syn_data$Y), exprs(median_normalization(ms)))
#   expect_equal(median_normalization(syn_data$Y), assay(median_normalization(as(ms, "SummarizedExperiment"))))
#
#
#   fit <- proDA(ms, ~ group, verbose=TRUE)
#
# })




test_that("subsampling works", {

  se <- generate_synthetic_data(n_proteins = 100, return_summarized_experiment = TRUE)

  fit <- proDA(se, ~ group, n_subsample = 50, verbose=TRUE)
  fit2 <- proDA(se, ~ group, verbose=TRUE)

  expect_gt(cor(c(predict(fit)), c(predict(fit2))), 0.95)

  # Check that there is no problem if n_subsample > nrow(se)
  fit3 <- proDA(se, ~ group, n_subsample = 200, verbose=TRUE)
})




