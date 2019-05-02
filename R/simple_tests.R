

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



#'
#'
#'
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




