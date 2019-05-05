
#' Predict the parameters or values of additional proteins
#'
#' This function can either predict the abundance matrix for proteins
#' (\code{type = "response"}) without missing values according to the
#' linear probabilistic dropout model, fitted with \code{proDA()}. Or, it
#' can predict the feature parameters for additional proteins given their
#' abundances including missing values after estimating the hyper-
#' parameters on a dataset with the same sample structure
#' (\code{type = "feature_parameters"}).
#'
#' \strong{Note:} this method behaves a little different from what one might
#' expect from the classical \code{predict.lm()} function, because
#' \code{object} is not just a single set of coefficients for one fit, but
#' many fits (ie. one for each protein) with some more hyper-parameters. The
#' classical \code{predict} function predicts the response for new samples.
#' This function does not support this, instead it is useful for getting a
#' matrix without missing values for additional proteins.
#'
#' @param object an 'proDAFit' object that is produced by \code{proDA()}.
#' @param newdata a matrix or a SummarizedExperiment which contains
#'   the new abundances for which values are predicted.
#' @param type either "response" or "feature_parameters". Default:
#'   \code{"response"}
#' @param ... additional parameters for the construction of the
#'   'proDAFit' object.
#'
#'
#' @return If \code{type = "response"} a matrix with the same dimensions
#'   as \code{object}. Or, if \code{type = "feature_parameters"} a
#'   'proDAFit' object with the same hyper-parameters and column data
#'   as \code{object}, but new fitted \code{rowData()}.
#'
#' @export
setMethod("predict", signature = "proDAFit", function(object, newdata,
                                                      type=c("response", "feature_parameters"), ...){
  type <- match.arg(type, c("response", "feature_parameters"))

  if(type == "response" && (missing(newdata) || is.null(newdata))){
    return(object$coefficients %*% t(object$design))
  }else if(type == "feature_parameters" && (missing(newdata) || is.null(newdata))){
    return(object)
  }

  # Check if newdata is applicable
  stopifnot(ncol(object) == ncol(newdata))
  stopifnot(colnames(object) == colnames(newdata))

  feat_data <- fit_feature_parameters(object, newdata)
  feat_df <- feat_data$feature_df
  coef_mat <- feat_data$coefficient_matrix

  new_fit_object <- proDAFit(newdata, NULL,
                      dropout_curve_position = object$hyper_parameters$dropout_curve_position,
                      dropout_curve_scale = object$hyper_parameters$dropout_curve_scale,
                      feature_parameters = feat_df,
                      coefficients = coef_mat,
                      design_matrix = design(object),
                      design_formula = design(object, formula=TRUE),
                      reference_level = reference_level(object),
                      location_prior_mean = object$hyper_parameters$location_prior_mean,
                      location_prior_scale = object$hyper_parameters$location_prior_scale,
                      location_prior_df = object$hyper_parameters$location_prior_df,
                      variance_prior_scale = object$hyper_parameters$variance_prior_scale,
                      variance_prior_df = object$hyper_parameters$variance_prior_df,
                      convergence = object$convergence, ...)


  if(type == "response"){
    predict(new_fit_object, type="response")
  }else{
    new_fit_object
  }
})


fit_feature_parameters <- function(fit, newdata){

  res_unreg <- lapply(seq_len(nrow(newdata)), function(i){
    pd_lm.fit(newdata[i, ], design(fit),
              dropout_curve_position = fit$hyper_parameters$dropout_curve_position,
              dropout_curve_scale = fit$hyper_parameters$dropout_curve_position)
  })
  if(! is.null(fit$hyper_parameters$location_prior_mean) ||
     ! is.null(fit$hyper_parameters$variance_prior_scale)){
    res_reg <- lapply(seq_len(nrow(newdata)), function(i){
      pd_lm.fit(newdata[i, ], design(fit),
                dropout_curve_position = fit$hyper_parameters$dropout_curve_position,
                dropout_curve_scale = fit$hyper_parameters$dropout_curve_position,
                location_prior_mean = fit$hyper_parameters$location_prior_mean,
                location_prior_scale = fit$hyper_parameters$location_prior_scale,
                variance_prior_scale = fit$hyper_parameters$variance_prior_scale,
                variance_prior_df = fit$hyper_parameters$variance_prior_df,
                location_prior_df = fit$hyper_parameters$location_prior_df)
    })
  }else{
    res_reg <- res_unreg
  }

  feat_df <- as.data.frame(mply_dbl(res_reg, function(f){
    unlist(f[-1])
  }, ncol = 4))
  # feat_df$rss <- vapply(res_unreg, function(x) x[["rss"]], 0.0)
  coef_mat <- mply_dbl(res_reg, function(f){
    f$coefficients
  }, ncol=ncol(design(fit)))
  colnames(coef_mat) <- names(res_reg[[1]]$coefficients)

  list(feature_df = feat_df,
       coefficient_matrix = coef_mat)
}
