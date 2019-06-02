


#' Get the abundance matrix
#'
#' @param object the object to get from
#' @param ... additional arguments used by the concrete implementation
#'
#' @examples
#'   syn_data <- generate_synthetic_data(n_proteins = 10)
#'   fit <- proDA(syn_data$Y, design = syn_data$groups)
#'   abundances(fit)
#'
#' @return the original matrix that was fitted
#'
#' @seealso \link{accessor_methods} for the implementation for a 'proDAFit' object
#' @export
setGeneric("abundances", function(object, ...) standardGeneric("abundances"))

#' Get the hyper parameters
#'
#' @param object the object to get from
#' @param ... additional arguments used by the concrete implementation
#'
#' @return a list with the values for each fitted hyper-parameter
#'
#' @examples
#'   syn_data <- generate_synthetic_data(n_proteins = 10)
#'   fit <- proDA(syn_data$Y, design = syn_data$groups)
#'   hyper_parameters(fit)
#'
#' @seealso \link{accessor_methods} for the implementation for a 'proDAFit' object
#' @export
setGeneric("hyper_parameters", function(object, ...) standardGeneric("hyper_parameters"))

#' Get the feature parameters
#'
#' @param object the object to get from
#' @param ... additional arguments used by the concrete implementation
#'
#' @return a data.frame with information on each fit
#'
#' @examples
#'   syn_data <- generate_synthetic_data(n_proteins = 10)
#'   fit <- proDA(syn_data$Y, design = syn_data$groups)
#'   feature_parameters(fit)
#'
#' @seealso \link{accessor_methods} for the implementation for a 'proDAFit' object
#' @export
setGeneric("feature_parameters", function(object, ...) standardGeneric("feature_parameters"))

#' Get the coefficients
#'
#' @param object the object to get from
#' @param ... additional arguments used by the concrete implementation
#'
#' @return a numeric matrix of size `nrow(fit) * p`
#'
#' @examples
#'   syn_data <- generate_synthetic_data(n_proteins = 10)
#'   fit <- proDA(syn_data$Y, design = syn_data$groups)
#'   coefficients(fit)
#'
#' @seealso \link{accessor_methods} for the implementation for a 'proDAFit' object
#' @export
setGeneric("coefficients", function(object, ...) standardGeneric("coefficients"))

#' Get the coefficients
#'
#' @param object the object to get from
#' @param ... additional arguments used by the concrete implementation
#'
#' @return a list with as many entries as rows in the data. Each element is
#'   a p*p matrix
#'
#' @examples
#'   syn_data <- generate_synthetic_data(n_proteins = 10)
#'   fit <- proDA(syn_data$Y, design = syn_data$groups)
#'   coefficient_variance_matrices(fit)
#'
#' @seealso \link{accessor_methods} for the implementation for a 'proDAFit' object
#' @export
setGeneric("coefficient_variance_matrices", function(object, ...) standardGeneric("coefficient_variance_matrices"))


#' Get the reference level
#'
#' @param object the object to get from
#' @param ... additional arguments used by the concrete implementation
#'
#' @return a string
#'
#' @examples
#'   syn_data <- generate_synthetic_data(n_proteins = 10)
#'   fit <- proDA(syn_data$Y, design = syn_data$groups)
#'   reference_level(fit)
#'
#' @seealso \link{accessor_methods} for the implementation for a 'proDAFit' object
#' @export
setGeneric("reference_level", function(object, ...) standardGeneric("reference_level"))


#' Get the convergence information
#'
#' @param object the object to get from
#' @param ... additional arguments used by the concrete implementation
#'
#' @return a list with information on the convergence
#'
#' @examples
#'   syn_data <- generate_synthetic_data(n_proteins = 10)
#'   fit <- proDA(syn_data$Y, design = syn_data$groups)
#'   convergence(fit)
#'
#' @seealso \link{accessor_methods} for the implementation for a 'proDAFit' object
#' @export
setGeneric("convergence", function(object, ...) standardGeneric("convergence"))


#' Calculate an approximate distance for 'object'
#'
#' @param object the object for which the distance is approximated
#' @param ... additional arguments used by the concrete implementation
#'
#' @return a list with two elements: `mean` and `sd` both are formally
#'   of class "dist"
#'
#' @examples
#'   syn_data <- generate_synthetic_data(n_proteins = 10)
#'   fit <- proDA(syn_data$Y, design = syn_data$groups)
#'   dist_approx(fit)
#'
#' @seealso \code{\link[stats]{dist}} for the base R function and
#'   \code{\link[=dist_approx_impl]{dist_approx()}} concrete implementation
#'   for 'proDAFit' objects
#' @export
setGeneric("dist_approx", function(object, ...) standardGeneric("dist_approx"))


#' Get the result_names
#'
#' @param fit the fit to get the result_names from
#' @param ... additional arguments used by the concrete implementation
#'
#' @return a character vector
#'
#' @examples
#'   syn_data <- generate_synthetic_data(n_proteins = 10)
#'   fit <- proDA(syn_data$Y, design = syn_data$groups)
#'   result_names(fit)
#'
#' @export
setGeneric("result_names", function(fit, ...) standardGeneric("result_names"))

