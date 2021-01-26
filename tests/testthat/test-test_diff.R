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


test_that("parse_contrast can handle problematic input", {
  
  expect_equal(parse_contrast(A, levels = c(LETTERS[1:5], "A:B")),
               c(A=1, B=0, C=0, D=0, E=0, `A:B`=0))
  
  expect_equal(parse_contrast(`A:B`, levels = c(LETTERS[1:5], "A:B")),
               c(A=0, B=0, C=0, D=0, E=0, `A:B`=1))
  
  expect_equal(parse_contrast(`A:B` - C, levels = c(LETTERS[1:5], "A:B")),
               c(A=0, B=0, C=-1, D=0, E=0, `A:B`=1))
  
  expect_equal(parse_contrast("`A:B`", levels = c(LETTERS[1:5], "A:B")),
               c(A=0, B=0, C=0, D=0, E=0, `A:B`=1))
  
  expect_equal(parse_contrast("`A:B` - C", levels = c(LETTERS[1:5], "A:B")),
               c(A=0, B=0, C=-1, D=0, E=0, `A:B`=1))
  
})

test_that("Parser contrast can handle reference to object in environment", {
  c1 <- parse_contrast(A - B, levels = LETTERS[1:2])
  c2 <- parse_contrast("A - B", levels = LETTERS[1:2])
  string <- "A - B"
  c3 <- parse_contrast(string, levels = LETTERS[1:2])
  c4 <- parse_contrast(paste0("A", " - ", "B"), levels = LETTERS[1:2])
  expect_equal(c1, c2)
  expect_equal(c1, c3)
  expect_equal(c1, c4)
})


test_that("parse contrast can handle being inside a function", {

  fnc <- function(){
    string <- "A - B"
    parse_contrast(string, levels = LETTERS[1:2])
  }

  c1 <- parse_contrast(A - B, levels = LETTERS[1:2])
  c2 <- fnc()
  expect_equal(c1, c2)
})


test_that("result names regex is approximately good", {
  regex <- "^[_.]?[[:alpha:]]+[[:alnum:]_.]*$"
  expect_true(grepl(regex, "tmp"))
  expect_true(grepl(regex, "t.mp"))
  expect_true(grepl(regex, "_tmp"))
  expect_true(grepl(regex, "y"))
  expect_true(grepl(regex, "s"))
  expect_true(grepl(regex, "A034245"))
  expect_true(grepl(regex, ".a00"))
  expect_false(grepl(regex, "."))
  expect_false(grepl(regex, "_"))
  expect_false(grepl(regex, "."))
  expect_false(grepl(regex, ""))
  expect_false(grepl(regex, ".0234"))
  expect_false(grepl(regex, "_42a"))
  expect_false(grepl(regex, "afdsa:fdsa"))
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
  set.seed(2)
  data <- matrix(rnorm(50 * 5, mean=20), ncol=5, nrow=50)
  data[invprobit(data, 19, -1) > runif(5 * 50)] <- NA
  colnames(data) <- paste0("sample_", LETTERS[seq_len(ncol(data))])
  rownames(data) <- paste0("protein", seq_len(nrow(data)))
  annot_df <- data.frame(Name = paste0("sample_", LETTERS[seq_len(ncol(data))]),
                         cond = c(rep(c("A", "B"), each=2), "B"),
                         num = runif(ncol(data)))
  se <- SummarizedExperiment(data, colData = annot_df)

  fit <- proDA(se, ~ cond + num, moderate_location = TRUE, moderate_variance = TRUE)
  test_res_f <- test_diff(fit, reduced_model = ~ num + 1, n_max = 1000, verbose=TRUE)
  test_res_wald <- test_diff(fit, condB)
  expect_gt(cor(test_res_wald$pval, test_res_f$pval, use="complete.obs"), 0.95)
})





test_that("pd_row_f_test works", {

  set.seed(1)
  se <- generate_synthetic_data(n_proteins = 10, n_conditions = 3,
                                n_replicates = 2,
                                return_summarized_experiment = TRUE)

  set.seed(1)
  fit <- proDA(se, ~ group + 1, max_iter = 10)
  std_test_res <- test_diff(fit, reduced_model = ~1)

  set.seed(1)
  f_test_res <- pd_row_f_test(se[ ,1:2,drop=FALSE],
                              se[ ,3:4,drop=FALSE],
                              se[ ,5:6,drop=FALSE],
                              max_iter = 10,
                              return_fit = TRUE)

  expect_equal(unlist(fit$hyper_parameters), unlist(f_test_res$fit$hyper_parameters),
               tolerance = 1e-3)
  expect_gt(cor(std_test_res$pval, f_test_res$test_results$pval, use = "complete.obs"), 0.9999)
  expect_equal(fit$coefficients[,1], f_test_res$fit$coefficients[,1], tolerance = 1e-3)
  expect_gt(cor(log(fit$feature_parameters$s2), log(f_test_res$fit$feature_parameters$s2)), 0.999)
  expect_equal(fit$coefficients%*% t(design(fit)), f_test_res$fit$coefficients%*% t(design(f_test_res$fit)), tolerance = 1e-3)
  Pred_var_wi <- mply_dbl(seq_len(nrow(se)), function(i){
    sapply(seq_len(ncol(se)), function(j)
      t(design(fit)[j,]) %*% fit$coefficient_variance_matrices[[i]] %*% design(fit)[j,])
  }, ncol=ncol(se))
  Pred_var_wio <- mply_dbl(seq_len(nrow(se)), function(i){
    sapply(seq_len(ncol(se)), function(j)
      t(design(f_test_res$fit)[j,]) %*% f_test_res$fit$coefficient_variance_matrices[[i]] %*% design(f_test_res$fit)[j,])
  }, ncol=ncol(se))
  expect_gt(cor(c(Pred_var_wi), c(Pred_var_wio)), 0.99)

})





test_that("Wald test works as good limma", {

  se <- generate_synthetic_data(n_proteins = 100, n_conditions = 4, dropout_curve_position = -100,
                                return_summarized_experiment = TRUE)

  fit <- proDA(se, ~ group, moderate_location = FALSE)
  test_res <- test_diff(fit, "groupCondition_3")

  dm <- design(fit)
  lim_fit <- limma::lmFit(assay(se,1), design = dm)
  lim_fit <- limma::contrasts.fit(lim_fit, limma::makeContrasts(groupCondition_3, levels = result_names(fit)))
  lim_fit <- limma::eBayes(lim_fit)

  test_lim <- limma::topTable(lim_fit, sort.by = "none", number = 100)
  expect_gt(cor(test_res$avg_abundance, test_lim$AveExpr), 0.99)
  expect_gt(cor(test_res$diff, test_lim$logFC), 0.99)
  expect_gt(cor(test_res$pval, test_lim$P.Value), 0.99)

})


test_that("Wald test works with missing values", {
  set.seed(2)
  se <- generate_synthetic_data(n_proteins = 1000, n_conditions = 4, dropout_curve_position = 17,
                                return_summarized_experiment = TRUE)

  fit <- proDA(se, ~ group, verbose = TRUE)
  test_res <- test_diff(fit, "groupCondition_3")

  dm <- design(fit)
  lim_fit <- limma::lmFit(assay(se,"full_observations"), design = dm)
  lim_fit <- limma::contrasts.fit(lim_fit, limma::makeContrasts(groupCondition_3, levels = result_names(fit)))
  lim_fit <- limma::eBayes(lim_fit)

  test_lim <- limma::topTable(lim_fit, sort.by = "none", number = 1000)
  expect_gt(cor(test_res$avg_abundance, test_lim$AveExpr), 0.9)
  expect_gt(cor(test_res$diff, test_lim$logFC), 0.7)
  expect_gt(cor(test_res$pval, test_lim$P.Value), 0.7)

})






