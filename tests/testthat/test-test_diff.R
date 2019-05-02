context("test-test_diff")

test_that("parse_contrast works", {

  expect_equal(parse_contrast(A, levels = LETTERS[1:5]),
               c(A=1, B=0, C=0, D=0, E=0))

  expect_equal(parse_contrast(A - B, levels = LETTERS[1:5]),
               c(A=1, B=-1, C=0, D=0, E=0))

  expect_equal(parse_contrast("A - B", levels = LETTERS[1:5]),
               c(A=1, B=-1, C=0, D=0, E=0))

  expect_equal(parse_contrast(A * 3, levels = LETTERS[1:5]),
               c(A=3, B=0, C=0, D=0, E=0))

  expect_equal(parse_contrast(B, levels = c("Intercept", LETTERS[2:5]), reference_level = "A"),
               c(Intercept = 0, B=1, C=0, D=0, E=0))

  expect_error(parse_contrast(B - A, levels = c("Intercept", LETTERS[2:5]), reference_level = "A"))

  expect_equal(parse_contrast("Intercept", levels = c("Intercept")),
               c(Intercept = 1))

})




test_that("F works", {

  y <- rnorm(10, mean=20)
  df <- data.frame(cond = c(rep("A", 6), rep("B", 4)),
                   num = rnorm(10))
  data <- matrix(y, nrow=1)
  colnames(data) <- paste0("sample_", LETTERS[seq_len(ncol(data))])
  rownames(data) <- paste0("protein", seq_len(nrow(data)))
  fit <- proDA(data, ~ cond + num, col_data=df,
               moderate_location = FALSE, moderate_variance = FALSE)
  glin_m <- glm(y ~ cond + num, data=df)
  # anova(glin_m, glm(y ~ 1), test="F")
  # test_diff(fit, reduced_model = ~ 1)
  expect_equal(anova(glin_m, glm(y ~ 1), test="F")[2, 6],  test_diff(fit, reduced_model = ~ 1)$pval)

})



test_that("F test works with missing data", {

  data <- matrix(rnorm(50 * 5, mean=20), ncol=5, nrow=50)
  data[invprobit(data, 19, -1) > runif(5 * 50)] <- NA
  colnames(data) <- paste0("sample_", LETTERS[seq_len(ncol(data))])
  rownames(data) <- paste0("protein", seq_len(nrow(data)))
  annot_df <- data.frame(Name = paste0("sample_", LETTERS[seq_len(ncol(data))]),
                         cond = c(rep(c("A", "B"), each=2), "B"),
                         num = runif(ncol(data)))
  se <- SummarizedExperiment(data, colData = annot_df)

  fit <- proDA(se, ~ cond + num, moderate_location = FALSE, moderate_variance = FALSE)
  test_diff(fit, condB, n_max = 10, verbose=TRUE)
  test_res <- test_diff(fit, reduced_model = ~ num, n_max = 1000, verbose=TRUE)
  head(test_res, n=10)

  expect_gt(cor(test_diff(fit, condB)$pval, test_res$pval, use="complete.obs"), 0.95)
})





test_that("pd_row_f_test works", {

  set.seed(2)
  n_proteins <- 100
  data <- matrix(rnorm(n_proteins * 6, mean=20), ncol=6, nrow=n_proteins)
  data[invprobit(data, 19, -1) > runif(prod(dim(data)))] <- NA
  colnames(data) <- paste0("sample_", LETTERS[seq_len(ncol(data))])
  rownames(data) <- paste0("protein", seq_len(nrow(data)))
  annot_df <- data.frame(Name = paste0("sample_", LETTERS[seq_len(ncol(data))]),
                         cond = c(rep(c("A", "B", "C"), each=2)),
                         num = runif(ncol(data)))
  se <- SummarizedExperiment(data, colData = annot_df)

  set.seed(1)
  fit <- proDA(se, ~ cond + 1, max_iter = 10,
               verbose = TRUE)
  std_test_res <- test_diff(fit, reduced_model = ~1)

  set.seed(1)
  f_test_res <- pd_row_f_test(data[ ,1:2,drop=FALSE],
                              data[ ,3:4,drop=FALSE],
                              data[ ,5:6,drop=FALSE],
                              max_iter = 10,
                              return_fit = TRUE, verbose=TRUE)

  expect_equal(unlist(fit$hyper_parameters), unlist(f_test_res$fit$hyper_parameters))
  # ggplot2::qplot(std_test_res$pval, f_test_res$test_results$pval,
  #                color= as.factor(f_test_res$test_results$n_obs))
  expect_gt(cor(std_test_res$pval, f_test_res$test_results$pval, use = "complete.obs"), 0.9999)

})










