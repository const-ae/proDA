

#' proDA Class Definition
#'
#' @export
.proDAFit <- setClass("proDAFit",
  slots = c(
    location_prior_mean = "numeric",
    location_prior_scale = "numeric",
    location_prior_df = "numeric",
    variance_prior_scale = "numeric",
    variance_prior_df = "numeric",
    design_matrix = "matrix",
    design_formula = "ANY",
    reference_level = "ANY",
    convergence = "list"
  ),
  contains = "SummarizedExperiment"
)

proDAFit <- function(data, col_data,
                     dropout_curve_position, dropout_curve_scale,
                     feature_parameters,
                     coefficients,
                     design_matrix, design_formula, reference_level,
                     location_prior_mean, location_prior_scale, location_prior_df,
                     variance_prior_scale, variance_prior_df,
                     convergence,
                     ...){

  se <- NULL
  if(is(data, "SummarizedExperiment")){
    assayNames(data)[1] <- "abundances"
    se <- data
    if(! is.null(col_data)){
      colData(se) <- cbind(colData(se), col_data)
    }
  }else if(is.matrix(data)){
    se <- SummarizedExperiment(assays=list(abundances=data), colData = col_data, ...)
  }else{
    stop("data must be a matrix or a SummarizedExperiment")
  }

  if(! is.numeric(dropout_curve_position) ||  length(dropout_curve_position) != ncol(se)){
    stop("dropout_curve_position must be numeric vector with one entry for each column")
  }
  if(! is.numeric(dropout_curve_scale) ||  length(dropout_curve_scale) != ncol(se)){
    stop("dropout_curve_position must be numeric vector with one entry for each column")
  }
  dropout_df <- DataFrame(dropout_curve_position, dropout_curve_scale)
  mcols(dropout_df) <- DataFrame(type = "hyper_parameter", description =
           c("The intensity where the chance to observe a protein is 50%",
             "How broad the sigmoidal dropout curve is."))
  colData(se) <- cbind(colData(se), dropout_df)


  if(! is.data.frame(feature_parameters) ||  !all(vapply(feature_parameters, is.numeric, FALSE)) ||
     nrow(feature_parameters) != nrow(se)){
    stop("feature_parameters must be a data.frame with as many rows as data and numeric columns")
  }
  feature_params_df <- DataFrame(feature_parameters)
  mcols(feature_params_df) <- DataFrame(type = "feature_parameter",
                                        description = "")
  rowData(se) <- cbind(rowData(se), feature_params_df)
  if(! is.matrix(coefficients) ||
     nrow(coefficients) != nrow(se)){
    stop("coefficients must be a martix with as many rows as data")
  }
  coefficients_df <- DataFrame(coefficients)
  mcols(coefficients_df) <- DataFrame(type = "coefficient",
                                        description = "The MAP estimate")
  rowData(se) <- cbind(rowData(se), coefficients_df)

  if(! is.matrix(design_matrix) || nrow(design_matrix) != ncol(se)){
    stop("design_matrix must be a matrix and the number of rows must match ",
         "the number of columns in data")
  }
  if(! is.null(design_formula) && ! inherits(design_formula, "formula")){
    stop("design_formula must be either NULL or of type formula")
  }
  if(! is.null(reference_level) && (! is.character(reference_level) || length(reference_level) != 1)){
    stop("reference_level must be a single string or NULL")
  }

  if(! is.numeric(location_prior_mean) || length(location_prior_mean) != 1){
    stop("location_prior_mean must be a single number")
  }
  if(! is.numeric(location_prior_scale) || length(location_prior_scale) != 1){
    stop("location_prior_scale must be a single number")
  }
  if(! is.numeric(location_prior_df) || length(location_prior_df) != 1){
    stop("location_prior_df must be a single number")
  }
  if(! is.numeric(variance_prior_scale) || length(variance_prior_scale) != 1){
    stop("variance_prior_scale must be a single number")
  }
  if(! is.numeric(variance_prior_df) || length(variance_prior_df) != 1){
    stop("variance_prior_df must be a single number")
  }

  if(any(names(convergence) != c("successful", "iterations", "error"))){
    stop("convergence must be a list with three elements: successful, iterations, error")
  }

  .proDAFit(se,
           location_prior_mean = location_prior_mean,
           location_prior_scale = location_prior_scale,
           location_prior_df = location_prior_df,
           variance_prior_scale = variance_prior_scale,
           variance_prior_df = variance_prior_df,
           design_matrix = design_matrix,
           design_formula = design_formula,
           reference_level = reference_level,
           convergence = convergence
           )


}


S4Vectors::setValidity2("proDAFit", function(object){

})

