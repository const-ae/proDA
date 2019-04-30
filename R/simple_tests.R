













#'
#'
#'
#' @export
pd_row_t_test <- function(X, Y,
                      moderate_location = TRUE,
                      moderate_variance = TRUE,
                      alternative = c("two.sided", "greater", "less"),
                      pval_adjust_method = "BH",
                      location_prior_df = 3,
                      max_iter = 20,
                      epsilon = 1e-3){

  design <- c(rep("Condition1", ncol(X)), rep("Condition2", ncol(Y)))
  data <- cbind(X, Y)
  fit <- proDA(data, design, moderate_location = moderate_location,
               moderate_variance = moderate_variance, max_iter = max_iter,
               epsilon = epsilon)
  test_diff(fit, Condition1 - Condition2,
            alternative = alternative,
            pval_adjust_method = pval_adjust_method)
}

