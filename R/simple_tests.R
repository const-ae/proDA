

#' Row-wise tests of difference using the probabilistic dropout model
#'
#' This is a helper function that combines the call of \code{proDA()}
#' and \code{test_diff()}. If you need more flexibility use those
#' functions.
#'
#' The \code{pd_row_t_test} is not acutally doing a t-test, but rather
#' a Wald test. But, as the two are closely related and term t-test is ≈
#' more widely understood, we choose to use that name.
#'
#' @param X,Y,... the matrices for condition 1, 2 and so on. They must
#'   have the same number of rows.
#' @param groups a factor or character vector with that assignes the
#'   columns of \code{X} to different conditions. This parameter is
#'   only applicable for the F-test and must be specified if only a
#'   single matrix is provided.
#' @param return_fit boolean that signals that in addition to the
#'   data.frame with the hypothesis test results, the fit from
#'   \code{proDA()} is returned. Default: \code{FALSE}
#' @inheritParams proDA
#' @inheritParams test_diff
#'
#' @return
#'   If \code{return_fit == FALSE} a data.frame is returned with the content
#'   that is described in \code{\link{test_diff}}.
#'
#'   If \code{return_fit == TRUE} a list is returned with two elements:
#'   \code{fit} with a reference to the object returned from \code{proDA()}
#'   and a \code{test_result()} with the data.frame returned from
#'   \code{test_diff()}.
#'
#'
#'
#' @seealso \code{\link{proDA}} and \code{\link{test_diff}} for more
#'   flexible versions. The  function was inspired
#'   by the \code{\link[genefilter]{rowttests}} function in the genefilter
#'   package.
#'
#'
#' @examples
#'   data1 <- matrix(rnorm(10 * 3), nrow=10)
#'   data2 <- matrix(rnorm(10 * 4), nrow=10)
#'   data3 <- matrix(rnorm(10 * 2), nrow=10)
#'
#'   # Comparing two datasets
#'   pd_row_t_test(data1, data2)
#'
#'   # Comparing multiple datasets
#'   pd_row_f_test(data1, data2, data3)
#'
#'   # Alternative
#'   data_comb <- cbind(data1, data2, data3)
#'   pd_row_f_test(data_comb,
#'      groups = c(rep("A",3), rep("B", 4), rep("C", 2)))
#'
#'   # t.test, lm, pd_row_t_test, and pd_row_f_test are
#'   # approximately equivalent on fully observed data
#'   set.seed(1)
#'   x <- rnorm(5)
#'   y <- rnorm(5, mean=0.3)
#'
#'   t.test(x, y)
#'   summary(lm(c(x, y) ~ cond,
#'              data = data.frame(cond = c(rep("x", 5),
#'                                         rep("y", 5)))))$coef[2,]
#'   pd_row_t_test(matrix(x, nrow=1), matrix(y, nrow=1),
#'                 moderate_location = FALSE,
#'                 moderate_variance = FALSE)
#'   pd_row_f_test(matrix(x, nrow=1), matrix(y, nrow=1),
#'                 moderate_location = FALSE,
#'                 moderate_variance = FALSE)
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
                      epsilon = 1e-3,
                      return_fit = FALSE,
                      verbose=FALSE){

  design <- c(rep("Condition1", ncol(X)), rep("Condition2", ncol(Y)))
  data <- cbind(X, Y)
  fit <- proDA(data, design, moderate_location = moderate_location,
               moderate_variance = moderate_variance, max_iter = max_iter,
               epsilon = epsilon, verbose=verbose)
  test_res <- test_diff(fit, Condition1 - Condition2,
            alternative = alternative,
            pval_adjust_method = pval_adjust_method,
            verbose = verbose)

  if(return_fit){
    list(fit = fit, test_results = test_res)
  }else{
    test_res
  }

}



#' @rdname pd_row_t_test
#' @export
pd_row_f_test <- function(X, ..., groups = NULL,
                          moderate_location = TRUE,
                          moderate_variance = TRUE,
                          pval_adjust_method = "BH",
                          location_prior_df = 3,
                          max_iter = 20,
                          epsilon = 1e-3,
                          return_fit = FALSE,
                          verbose=FALSE){
  additional_matrices <- list(...)

  if(length(additional_matrices) == 0 && (is.null(groups) ||
                            length(levels(as.factor(groups))) == 1) ){
    stop("Data for only one condition was specified. Please provide either a multiple ",
         "data matrices (pd_row_f_test(X, Y, Z, ...)) or provide the groups parameter",
         "(pd_row_f_test(X, groups=c('a', 'a', 'b', 'b', ...)).")
  }
  if(length(additional_matrices) > 0 && ! is.null(groups)){
    stop("Please specify ... or annotate the samples in X using groups, but not both")
  }
  if(length(additional_matrices) > 0){
    groups <- c(rep(paste0("condition_", 1), ncol(X)),
                unlist(lapply(seq_along(additional_matrices), function(idx){
      rep(paste0("condition_", idx+1), ncol(additional_matrices[[idx]]))
    })))
    X <- cbind(X, do.call(cbind, additional_matrices))
  }
  if(length(groups) != ncol(X)){
    stop("Length of groups must match the number of columns in X")
  }

  fit <- proDA(X, design = groups, moderate_location = moderate_location,
               moderate_variance = moderate_variance, max_iter = max_iter,
               epsilon = epsilon,verbose = verbose)


  test_res <- test_diff(fit, reduced_model = ~ 1, pval_adjust_method = pval_adjust_method,
            verbose = verbose)

  if(return_fit){
    list(fit = fit, test_results = test_res)
  }else{
    test_res
  }
}




