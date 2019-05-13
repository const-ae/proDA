context("test-proDAFit_S4")


test_that("proDAFit object construction works", {
  data <- matrix(rnorm(100), nrow=20, ncol=5)
  data[sample(1:100, 20)] <- NA
  colnames(data) <- paste0("Sample", 1:5)
  rownames(data) <- paste0("Protein_", 1:20)
  annot_data <- data.frame(sample = paste0("sample_", 1:5),
                           condition = c("A", "A", "B", "B", "C"))
  row_annot_data <- data.frame(name = paste0("protein_", 1:20))
  se <- SummarizedExperiment(assays = list(data=data), colData = annot_data, rowData = row_annot_data)

  rho <- rnorm(5, mean=20)
  zeta <- rnorm(5, mean=-1, sd=0.1)

  feat_params <- data.frame(n_approx = rgamma(20, 4),
                            df = 2,
                            s2 = extraDistr::rinvchisq(20, nu=2, tau=0.03),
                            n_obs = rpois(20, 8)
                            )
  coef_mat <- matrix(rnorm(40), nrow=20, ncol=2)
  colnames(coef_mat) <- c("Intercept", "beta1")
  coef_var <- lapply(seq_len(20), function(idx){
    diag(2)
  })


  conv <- list(successful = TRUE, iterations = 7, error =1e-5)

  dm <- stats::model.matrix.default(~ annot_data$condition)
  pf <- proDAFit(se,
           col_data = NULL,
           dropout_curve_position = rho,
           dropout_curve_scale = zeta,
           feature_parameters = feat_params,
           coefficients = coef_mat,
           coef_var = coef_var,
           design_matrix = dm,
           design_formula = NULL,
           reference_level = NULL,
           location_prior_mean = 20,
           location_prior_scale = 4.5,
           location_prior_df = 3,
           variance_prior_scale = 0.05,
           variance_prior_df = 2.1,
           convergence = conv
  )

  pf2 <- proDAFit(data,
                  col_data = annot_data,
                  dropout_curve_position = rho,
                  dropout_curve_scale = zeta,
                  feature_parameters = feat_params,
                  coefficients = coef_mat,
                  coef_var = coef_var,
                  design_matrix = dm,
                  design_formula = NULL,
                  reference_level = NULL,
                  location_prior_mean = 20,
                  location_prior_scale = 4.5,
                  location_prior_df = 3,
                  variance_prior_scale = 0.05,
                  variance_prior_df = 2.1,
                  convergence = conv
  )

  rowData(pf2) <- cbind(row_annot_data, rowData(pf2))
  expect_equal(pf, pf2)

  expect_equal(hyper_parameters(pf), pf$hyper_parameters)
  expect_equal(colData(pf), pf$colData)
  expect_error(pf$nrow)
})


